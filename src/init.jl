value(x) = x
cuify(x) = error("To use LinSolveGPUFactorize, you must do `using CuArrays`")
promote_u0(u0,p,t0) = u0
promote_tspan(u0,p,tspan,prob,kwargs) = tspan
get_tmp(x) = nothing
isdistribution(u0) = false

function SciMLBase.tmap(args...)
  error("Zygote must be added to differentiate Zygote? If you see this error, report it.")
end

const IS_OPENBLAS = Ref(true)

function __init__()
  @static if VERSION < v"1.7beta"
    blas = BLAS.vendor()
    IS_OPENBLAS[] = blas == :openblas64 || blas == :openblas
  else
    IS_OPENBLAS[] = occursin("openblas", BLAS.get_config().loaded_libs[1].libname)
  end
  @require ApproxFun="28f2ccd6-bb30-5033-b560-165f7b14dc2f" begin
    eval_u0(u0::ApproxFun.Fun) = false
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
    value(x::MonteCarloMeasurements.AbstractParticles) = mean(x.particles)

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
    include("tracker.jl")
  end

  @require CUDA="052768ef-5323-5732-b1bb-66c8b64840ba" begin
    cuify(x::AbstractArray) = CUDA.CuArray(x)
    default_factorize(A::CUDA.CuArray) = qr(A)

    ODE_DEFAULT_NORM(u::CUDA.CuArray,t) = sqrt(real(sum(abs2,u))/length(u))

    @inline function ODE_DEFAULT_NORM(u::CUDA.CuArray{<:ForwardDiff.Dual},t)
      sqrt(sum(abs2,value.(u)) / length(u))
    end

    @inline function ODE_DEFAULT_NORM(u::CUDA.CuArray{<:ForwardDiff.Dual},t::ForwardDiff.Dual)
      sqrt(sum(abs2,u) / length(u))
    end
  end

  @require ReverseDiff="37e2e3b7-166d-5795-8a7a-e32c996b4267" begin
    include("reversediff.jl")
  end

  @require Zygote="e88e6eb3-aa80-5325-afca-941959d7151f" begin
    function ∇tmap(cx, f, args...)
      ys_and_backs = SciMLBase.tmap((args...) -> Zygote._pullback(cx, f, args...), args...)
      if isempty(ys_and_backs)
        ys_and_backs, _ -> (NoTangent(),NoTangent())
      else
        ys, backs = Zygote.unzip(ys_and_backs)
        function ∇tmap_internal(Δ)
          Δf_and_args_zipped = SciMLBase.tmap((f, δ) -> f(δ), backs, Δ)
          Δf_and_args = Zygote.unzip(Δf_and_args_zipped)
          Δf = reduce(Zygote.accum, Δf_and_args[1])
          (Δf, Δf_and_args[2:end]...)
        end
        ys,∇tmap_internal
      end
    end

    function ∇responsible_map(cx, f, args...)
      ys_and_backs = SciMLBase.responsible_map((args...) -> Zygote._pullback(cx, f, args...), args...)
      if isempty(ys_and_backs)
        ys_and_backs, _ -> (NoTangent(),NoTangent())
      else
        ys, backs = Zygote.unzip(ys_and_backs)
        ys, function ∇responsible_map_internal(Δ)
          # Apply pullbacks in reverse order. Needed for correctness if `f` is stateful.
          Δf_and_args_zipped = SciMLBase.responsible_map((f, δ) -> f(δ), Zygote._tryreverse(SciMLBase.responsible_map, backs, Δ)...)
          Δf_and_args = Zygote.unzip(Zygote._tryreverse(SciMLBase.responsible_map, Δf_and_args_zipped))
          Δf = reduce(Zygote.accum, Δf_and_args[1])
          (Δf, Δf_and_args[2:end]...)
        end
      end
    end
  end

  @require GeneralizedGenerated="6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb" begin
    SciMLBase.numargs(::GeneralizedGenerated.RuntimeFn{Args}) where Args = GeneralizedGenerated.from_type(Args) |> length
  end
end
