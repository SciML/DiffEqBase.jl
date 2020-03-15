value(x) = x
cuify(x) = error("To use LinSolveGPUFactorize, you must do `using CuArrays`")
promote_u0(u0,p,t0) = u0
promote_tspan(u0,p,tspan,prob,kwargs) = tspan

if VERSION < v"1.4.0-DEV.635"
  # Piracy, should get upstreamed
  LinearAlgebra.ldiv!(Y::AbstractArray, A::AbstractArray, B::AbstractArray) = (copyto!(Y,B); ldiv!(A,Y))
  LinearAlgebra.ldiv!(Y::AbstractArray, A::AbstractMatrix, B::AbstractArray) = (copyto!(Y,B); ldiv!(A,Y))
  LinearAlgebra.ldiv!(x::AbstractArray,A::Diagonal,b::AbstractArray) = (x .= A.diag .\ b)
end

function __init__()
  @require ApproxFun="28f2ccd6-bb30-5033-b560-165f7b14dc2f" begin
    eval_u0(u0::ApproxFun.Fun) = false
  end

  @require Distributions="31c24e10-a181-5473-b8eb-7969acd0382f" begin
    handle_distribution_u0(_u0::Distributions.Sampleable) = rand(_u0)
  end

  @require ForwardDiff="f6369f11-7733-5829-9624-2563aa707210" begin

    promote_u0(u0::AbstractArray{<:ForwardDiff.Dual},p::AbstractArray{<:ForwardDiff.Dual},t0) = u0
    promote_u0(u0,p::AbstractArray{<:ForwardDiff.Dual},t0) = eltype(p).(u0)

    function promote_tspan(u0::AbstractArray{<:ForwardDiff.Dual},p,tspan::Tuple{<:ForwardDiff.Dual,<:ForwardDiff.Dual},prob,kwargs)
      return tspan
    end

    function promote_tspan(u0::AbstractArray{<:ForwardDiff.Dual},p,tspan,prob,kwargs)
      if (haskey(kwargs,:callback) && has_continuous_callback(kwargs[:callback])) ||
         (haskey(prob.kwargs,:callback) && has_continuous_callback(prob.kwargs[:callback]))

        return eltype(u0).(tspan)
      else
        return tspan
      end
    end

    value(x::Type{ForwardDiff.Dual{T,V,N}}) where {T,V,N} = V
    value(x::ForwardDiff.Dual) = value(ForwardDiff.value(x))

    @inline fastpow(x::ForwardDiff.Dual, y::ForwardDiff.Dual) = x^y

    # Support adaptive with non-dual time
    @inline ODE_DEFAULT_NORM(u::AbstractArray{<:ForwardDiff.Dual},::Any) = sqrt(sum(UNITLESS_ABS2∘value,u) / length(u))
    @inline ODE_DEFAULT_NORM(u::ForwardDiff.Dual,::Any) = abs(value(u))

    # When time is dual, it shouldn't drop the duals for adaptivity
    @inline ODE_DEFAULT_NORM(u::AbstractArray{<:ForwardDiff.Dual},::ForwardDiff.Dual) = sqrt(sum(UNITLESS_ABS2,u) / length(u))
    @inline ODE_DEFAULT_NORM(u::ForwardDiff.Dual,::ForwardDiff.Dual) = abs(u)

    # Type piracy. Should upstream
    Base.nextfloat(d::ForwardDiff.Dual{T,V,N}) where {T,V,N} = ForwardDiff.Dual{T}(nextfloat(d.value), d.partials)
    Base.prevfloat(d::ForwardDiff.Dual{T,V,N}) where {T,V,N} = ForwardDiff.Dual{T}(prevfloat(d.value), d.partials)

    struct DiffCache{T<:AbstractArray, S<:AbstractArray}
        du::T
        dual_du::S
    end

    function DiffCache(u::AbstractArray{T}, siz, ::Type{Val{chunk_size}}) where {T, chunk_size}
        DiffCache(u, zeros(ForwardDiff.Dual{nothing,T,chunk_size}, siz...))
    end

    dualcache(u::AbstractArray, N=Val{ForwardDiff.pickchunksize(length(u))}) = DiffCache(u, size(u), N)

    get_tmp(dc::DiffCache, u::AbstractArray{T}) where T<:ForwardDiff.Dual = reinterpret(T, dc.dual_du)
    get_tmp(dc::DiffCache, u::AbstractArray) = dc.du
  end

  @require Measurements="eff96d63-e80a-5855-80a2-b1b0885c5ab7" begin

    promote_u0(u0::AbstractArray{<:Measurements.Measurement},p::AbstractArray{<:Measurements.Measurement},t0) = u0
    promote_u0(u0,p::AbstractArray{<:Measurements.Measurement},t0) = eltype(p).(u0)

    value(x::Type{Measurements.Measurement{T}}) where {T} = T
    value(x::Measurements.Measurement) = Measurements.value(x)

    @inline fastpow(x::Measurements.Measurement, y::Measurements.Measurement) = x^y

    # Support adaptive steps should be errorless
    @inline function ODE_DEFAULT_NORM(u::AbstractArray{<:Measurements.Measurement,N},t) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(t))) / length(u))
    end
    @inline function ODE_DEFAULT_NORM(u::Array{<:Measurements.Measurement,N},t) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(t))) / length(u))
    end
    @inline ODE_DEFAULT_NORM(u::Measurements.Measurement,t) = abs(Measurements.value(u))
  end

  @require MonteCarloMeasurements="0987c9cc-fe09-11e8-30f0-b96dd679fdca" begin

    promote_u0(u0::AbstractArray{<:MonteCarloMeasurements.AbstractParticles},p::AbstractArray{<:MonteCarloMeasurements.AbstractParticles},t0) = u0
    promote_u0(u0,p::AbstractArray{<:MonteCarloMeasurements.AbstractParticles},t0) = eltype(p).(u0)

    value(x::Type{MonteCarloMeasurements.AbstractParticles{T,N}}) where {T,N} = T
    value(x::MonteCarloMeasurements.AbstractParticles) = mean(x)

    @inline fastpow(x::MonteCarloMeasurements.AbstractParticles, y::MonteCarloMeasurements.AbstractParticles) = x^y

    # Support adaptive steps should be errorless
    @inline function ODE_DEFAULT_NORM(u::AbstractArray{<:MonteCarloMeasurements.AbstractParticles,N},t) where {N}
      sqrt(mean(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(t))))
    end
    @inline function ODE_DEFAULT_NORM(u::AbstractArray{<:MonteCarloMeasurements.AbstractParticles,N},t::AbstractArray{<:MonteCarloMeasurements.AbstractParticles,N}) where {N}
      sqrt(mean(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(value.(t)))))
    end
    @inline ODE_DEFAULT_NORM(u::MonteCarloMeasurements.AbstractParticles,t) = abs(value(u))
  end

  @require Unitful="1986cc42-f94f-5a68-af5c-568840ba703d" begin
    # Support adaptive errors should be errorless for exponentiation
    value(x::Type{Unitful.AbstractQuantity{T,D,U}}) where {T,D,U} = T
    value(x::Unitful.AbstractQuantity) = x.val
    @inline function ODE_DEFAULT_NORM(u::AbstractArray{<:Unitful.AbstractQuantity,N},t) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(t))) / length(u))
    end
    @inline function ODE_DEFAULT_NORM(u::Array{<:Unitful.AbstractQuantity,N},t) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(t))) / length(u))
    end
    @inline ODE_DEFAULT_NORM(u::Unitful.AbstractQuantity,t) = abs(value(u))
    @inline UNITLESS_ABS2(x::Unitful.AbstractQuantity) = real(abs2(x)/oneunit(x)*oneunit(x))
  end

  @require Tracker="9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c" begin
    value(x::Type{Tracker.TrackedReal{T}}) where T = T
    value(x::Type{Tracker.TrackedArray{T,N,A}}) where {T,N,A} = Array{T,N}
    value(x::Tracker.TrackedReal)  = x.data
    value(x::Tracker.TrackedArray) = x.data

    @inline fastpow(x::Tracker.TrackedReal, y::Tracker.TrackedReal) = x^y
    @inline Base.any(f::Function,x::Tracker.TrackedArray) = any(f,Tracker.data(x))

    # Support adaptive with non-tracked time
    @inline function ODE_DEFAULT_NORM(u::Tracker.TrackedArray,t) where {N}
      sqrt(sum(abs2,value(u)) / length(u))
    end
    @inline function ODE_DEFAULT_NORM(u::AbstractArray{<:Tracker.TrackedReal,N},t) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(t))) / length(u))
    end
    @inline function ODE_DEFAULT_NORM(u::Array{<:Tracker.TrackedReal,N},t) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(t))) / length(u))
    end
    @inline ODE_DEFAULT_NORM(u::Tracker.TrackedReal,t) = abs(value(u))

    # Support TrackedReal time, don't drop tracking on the adaptivity there
    @inline function ODE_DEFAULT_NORM(u::Tracker.TrackedArray,t::Tracker.TrackedReal) where {N}
      sqrt(sum(abs2,u) / length(u))
    end
    @inline function ODE_DEFAULT_NORM(u::AbstractArray{<:Tracker.TrackedReal,N},t::Tracker.TrackedReal) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip(u,Iterators.repeated(t))) / length(u))
    end
    @inline function ODE_DEFAULT_NORM(u::Array{<:Tracker.TrackedReal,N},t::Tracker.TrackedReal) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip(u,Iterators.repeated(t))) / length(u))
    end
    @inline ODE_DEFAULT_NORM(u::Tracker.TrackedReal,t::Tracker.TrackedReal) = abs(u)
  end

  # Piracy, should get upstreamed
  @require CuArrays="3a865a2d-5b23-5a0f-bc46-62713ec82fae" begin
    cuify(x::AbstractArray) = CuArrays.CuArray(x)
    function LinearAlgebra.ldiv!(x::CuArrays.CuArray,_qr::CuArrays.CUSOLVER.CuQR,b::CuArrays.CuArray)
      _x = UpperTriangular(_qr.R) \ (_qr.Q' * reshape(b,length(b),1))
      x .= vec(_x)
      CuArrays.unsafe_free!(_x)
    end
    function findall_events(affect!,affect_neg!,prev_sign::CuArrays.CuArray,next_sign::CuArrays.CuArray)
      hasaffect::Bool = affect! !== nothing
      hasaffectneg::Bool = affect_neg! !== nothing
      f = (p,n)-> ((p < 0 && hasaffect) || (p > 0 && hasaffectneg)) && p*n<=0
      A = map(f,prev_sign,next_sign)
      out = findall(A)
      CuArrays.unsafe_free!(A)
      out
    end

    ODE_DEFAULT_NORM(u::CuArrays.CuArray,t) = sqrt(real(sum(abs2,u))/length(u))

    @require ForwardDiff="f6369f11-7733-5829-9624-2563aa707210" begin
      @inline function ODE_DEFAULT_NORM(u::CuArrays.CuArray{<:ForwardDiff.Dual},t)
        sqrt(sum(abs2,value.(u)) / length(u))
      end

      @inline function ODE_DEFAULT_NORM(u::CuArrays.CuArray{<:ForwardDiff.Dual},t::ForwardDiff.Dual)
        sqrt(sum(abs2,u) / length(u))
      end
    end
  end

  @require GeneralizedGenerated="6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb" begin
    numargs(::GeneralizedGenerated.RuntimeFn{Args}) where Args = GeneralizedGenerated.from_type(Args) |> length
  end

  @require Pardiso="46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2" begin
    mutable struct MKLPardisoFactorize
      A
      ps::Pardiso.MKLPardisoSolver
      verbose::Bool
      firsttime::Bool
    end
    MKLPardisoFactorize(;verbose=false) = MKLPardisoFactorize(nothing,Pardiso.MKLPardisoSolver(),verbose,true)
    function (p::MKLPardisoFactorize)(x,A,b,update_matrix=false;kwargs...)
      if p.firsttime
        Pardiso.set_phase!(p.ps, Pardiso.ANALYSIS)
        Pardiso.pardiso(p.ps, x, A, b)
        p.firsttime = false
      end

      if update_matrix
        Pardiso.set_phase!(p.ps, Pardiso.NUM_FACT)
        Pardiso.pardiso(p.ps, x, A, b)
        p.A = A
      end

      Pardiso.set_phase!(p.ps, Pardiso.SOLVE_ITERATIVE_REFINE)
      Pardiso.pardiso(p.ps, x, A, b)
    end
    function (p::MKLPardisoFactorize)(::Type{Val{:init}},f,u0_prototype)
      if eltype(u0_prototype) <: Complex
        mattype = Pardiso.COMPLEX_NONSYM
      else
        mattype = Pardiso.REAL_NONSYM
      end
      Pardiso.set_matrixtype!(p.ps, mattype)
      if p.verbose
          Pardiso.set_msglvl!(p.ps, Pardiso.MESSAGE_LEVEL_ON)
      end
      MKLPardisoFactorize(nothing,p.ps,p.verbose,true)
    end

    mutable struct PardisoFactorize
      A
      ps::Pardiso.PardisoSolver
      verbose::Bool
      firsttime::Bool
    end
    PardisoFactorize(;verbose=false) = PardisoFactorize(nothing,Pardiso.PardisoSolver(),verbose,true)
    function (p::PardisoFactorize)(x,A,b,update_matrix=false;kwargs...)
      if p.firsttime
        Pardiso.set_phase!(p.ps, Pardiso.ANALYSIS)
        Pardiso.pardiso(p.ps, x, A, b)
        p.firsttime = false
      end

      if update_matrix
        Pardiso.set_phase!(p.ps, Pardiso.NUM_FACT)
        Pardiso.pardiso(p.ps, x, A, b)
        p.A = A
      end

      Pardiso.set_phase!(p.ps, Pardiso.SOLVE_ITERATIVE_REFINE)
      Pardiso.pardiso(p.ps, x, A, b)
    end
    function (p::PardisoFactorize)(::Type{Val{:init}},f,u0_prototype)
      if eltype(u0_prototype) <: Complex
        mattype = Pardiso.COMPLEX_NONSYM
      else
        mattype = Pardiso.REAL_NONSYM
      end
      Pardiso.set_matrixtype!(p.ps, mattype)
      if p.verbose
          Pardiso.set_msglvl!(p.ps, Pardiso.MESSAGE_LEVEL_ON)
      end
      PardisoFactorize(nothing,p.ps,p.verbose,true)
    end

    mutable struct PardisoIterate
      A
      ps::Pardiso.PardisoSolver
      verbose::Bool
      firsttime::Bool
    end
    PardisoIterate(;verbose=false) = PardisoIterate(nothing,Pardiso.PardisoSolver(),verbose,true)
    function (p::PardisoIterate)(x,A,b,update_matrix=false;kwargs...)
      if p.firsttime
        Pardiso.set_phase!(p.ps, Pardiso.ANALYSIS)
        Pardiso.pardiso(p.ps, x, A, b)
        p.firsttime = false
      end

      #=
      # Pardiso iterative solver doesn't need a separate factorization phase
      # It will auto-update as it needs to.
      if update_matrix
        Pardiso.set_phase!(p.ps, Pardiso.NUM_FACT)
        Pardiso.pardiso(p.ps, x, A, b)
        p.A = A
      end
      =#

      Pardiso.set_phase!(ps, Pardiso.NUM_FACT_SOLVE_REFINE)
      Pardiso.pardiso(ps, X, A, B)
    end
    function (p::PardisoIterate)(::Type{Val{:init}},f,u0_prototype)
      if eltype(u0_prototype) <: Complex
        mattype = Pardiso.COMPLEX_NONSYM
      else
        mattype = Pardiso.REAL_NONSYM
      end
      Pardiso.set_matrixtype!(p.ps, mattype)
      Pardiso.set_solver!(ps, Pardiso.ITERATIVE_SOLVER)
      if p.verbose
          Pardiso.set_msglvl!(p.ps, Pardiso.MESSAGE_LEVEL_ON)
      end
      PardisoIterate(nothing,p.ps,p.verbose,true)
    end
  end
end
