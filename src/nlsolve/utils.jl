"""
    qrdelete!(Q, R, k)

Delete the left-most column of F = Q[:, 1:k] * R[1:k, 1:k] by updating Q and R.
Only Q[:, 1:(k-1)] and R[1:(k-1), 1:(k-1)] are valid on exit.
"""
function qrdelete!(Q::AbstractMatrix, R::AbstractMatrix, k::Int)
  n, m = size(Q)
  m == LinearAlgebra.checksquare(R) || throw(DimensionMismatch())
  1 ≤ k ≤ m || throw(ArgumentError())

  # apply Givens rotations
  for i in 2:k
      g = first(givens(R, i - 1, i, i))
      lmul!(g, R)
      rmul!(Q, g')
  end

  # move columns of R
  @inbounds for j in 1:(k-1)
    for i in 1:(k-1)
      R[i, j] = R[i, j + 1]
    end
  end

  Q, R
end

"""
    qradd!(Q, R, v, k)

Replace the right-most column of F = Q[:, 1:k] * R[1:k, 1:k] with v by updating Q and R.
This implementation modifies vector v as well. Only Q[:, 1:k] and R[1:k, 1:k] are valid on
exit.
"""
function qradd!(Q::AbstractMatrix, R::AbstractMatrix, v::AbstractVector, k::Int)
  n, m = size(Q)
  n == length(v) || throw(DimensionMismatch())
  m == LinearAlgebra.checksquare(R) || throw(DimensionMismatch())
  1 ≤ k ≤ m || throw(ArgumentError())

  @inbounds for i in 1:(k-1)
    q = view(Q, :, i)
    r = dot(q, v)

    R[i, k] = r
    axpy!(-r, q, v)
  end

  @inbounds begin
    d = norm(v)
    R[k, k] = d
    @.. @view(Q[:, k]) = v / d
  end

  Q, R
end

function qradd!(Q::AbstractMatrix, R::AbstractMatrix, v::Number, k::Int)
  1 == LinearAlgebra.checksquare(Q) == LinearAlgebra.checksquare(R) ||
    throw(DimensionMismatch())
  k == 1 || throw(ArgumentError())

  R[1, 1] = abs(v)
  Q[1, 1] = one(v)

  Q, R
end

function nlsolve_f end
function iip_get_uf end
function oop_get_uf end
function build_jac_config end
function resize_jac_config! end

# No J version
function iipnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  iipnlsolve(alg,u,uprev,p,t,dt,f,W,nothing,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
end

function iipnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  # define additional fields of cache of non-linear solver
  z = similar(u); dz = similar(u); tmp = similar(u); b = similar(u)
  k = zero(rate_prototype)

  uTolType = real(uBottomEltypeNoUnits)

  # create cache of non-linear solver
  if alg.nlsolve isa NLNewton
    nf = nlsolve_f(f, alg)

    if islinear(f)
      du1 = rate_prototype
      uf = nothing
      jac_config = nothing
      linsolve = alg.linsolve(Val{:init},nf,u)
    else
      du1 = zero(rate_prototype)
      # if the algorithm specializes on split problems the use `nf`
      # we pass this `alg` here just for identification purpose, because get_uf would be overloaded in different repos
      uf = iip_get_uf(alg,nf,t,p)
      jac_config = build_jac_config(alg,nf,uf,du1,uprev,u,tmp,dz)
      linsolve = alg.linsolve(Val{:init},uf,u)
    end
    # TODO: check if the solver is iterative
    weight = similar(u)

    nlcache = DiffEqBase.NLNewtonCache(true,W,J,dt,du1,uf,jac_config,linsolve,weight)
  elseif alg.nlsolve isa NLFunctional
    z₊ = similar(z)

    nlcache = DiffEqBase.NLFunctionalCache(z₊)
  elseif alg.nlsolve isa NLAnderson
    z₊ = similar(z)

    max_history = min(alg.nlsolve.max_history, alg.nlsolve.max_iter, length(z))
    Δz₊s = [zero(z) for i in 1:max_history]
    Q = Matrix{uEltypeNoUnits}(undef, length(z), max_history)
    R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
    γs = Vector{uEltypeNoUnits}(undef, max_history)
    dzold = zero(z)
    z₊old = zero(z)

    nlcache = DiffEqBase.NLAndersonCache(z₊,dzold,z₊old,Δz₊s,Q,R,γs)
  end

  # create non-linear solver
  nlsolver = NLSolver{typeof(alg.nlsolve),true,typeof(u),typeof(k),uTolType,typeof(γ),typeof(c),typeof(nlcache)}(z,dz,tmp,b,k,alg.nlsolve,one(uTolType),γ,c,10000,Convergence,nlcache)
end

# No J version
function oopnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  oopnlsolve(alg,u,uprev,p,t,dt,f,W,nothing,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)  
end

function oopnlsolve(alg,u,uprev,p,t,dt,f,W,J,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  # define additional fields of cache of non-linear solver (all aliased)
  z = uprev; dz = z; tmp = z; b = z; k = rate_prototype

  uTolType = real(uBottomEltypeNoUnits)

  # create cache of non-linear solver
  if alg.nlsolve isa NLNewton
    nf = nlsolve_f(f, alg)
    # only use `nf` if the algorithm specializes on split eqs
    uf = oop_get_uf(alg,nf,t,p)
    nlcache = DiffEqBase.NLNewtonConstantCache(true,W,J,dt,uf)
  elseif alg.nlsolve isa NLFunctional
    nlcache = DiffEqBase.NLFunctionalConstantCache()
  elseif alg.nlsolve isa NLAnderson
    max_history = min(alg.nlsolve.max_history, alg.nlsolve.max_iter, length(z))
    Δz₊s = Vector{typeof(z)}(undef, max_history)
    Q = Matrix{uEltypeNoUnits}(undef, length(z), max_history)
    R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
    γs = Vector{uEltypeNoUnits}(undef, max_history)

    nlcache = DiffEqBase.NLAndersonConstantCache(Δz₊s,Q,R,γs)
  end

  # create non-linear solver
  nlsolver = NLSolver{typeof(alg.nlsolve),false,typeof(u),typeof(k),uTolType,typeof(γ),typeof(c),typeof(nlcache)}(z,dz,tmp,b,k,alg.nlsolve,one(uTolType),γ,c,10000,Convergence,nlcache)
end

## resize NLSolver

function nlsolve_resize!(integrator::DEIntegrator, i::Int)
  (isdefined(integrator.alg, :nlsolve) && isdefined(integrator.cache, :nlsolver)) || return

  nlsolver = integrator.cache.nlsolver
  if nlsolver isa AbstractArray
    for idx in eachindex(nlsolver)
      resize!(nlsolver[idx], i)
    end
  else
    resize!(nlsolver, i)
  end

  nothing
end

function Base.resize!(nlsolver::NLSolver, i::Int)
  @unpack z,dz,tmp,ztmp,k,cache = nlsolver

  resize!(z, i)
  resize!(dz, i)
  resize!(tmp, i)
  resize!(ztmp, i)
  resize!(k, i)

  if nlsolver.alg isa NLAnderson
    resize!(cache, nlsolver.alg, i)
  else
    resize!(cache, i)
  end
end

Base.resize!(::AbstractNLSolverCache, ::Int) = nothing

function Base.resize!(nlcache::NLFunctionalCache, i::Int)
  resize!(nlcache.z₊, i)
  nothing
end

function Base.resize!(nlcache::NLNewtonCache, i::Int)
  @unpack du1,jac_config,linsolve,weight = nlcache

  resize!(du1, i)
  if jac_config !== nothing
    nlsolver.jac_config = resize_jac_config!(jac_config, i)
  end
  resize!(weight, i)

  nothing
end

function Base.resize!(nlcache::NLAndersonCache, nlalg::NLAnderson, i::Int)
  @unpack z₊, dzold, z₊old, γs, Δz₊s = nlcache

  resize!(z₊, i)
  resize!(dzold, i)
  resize!(z₊old, i)

  # determine new maximum history
  max_history_old = length(Δz₊s)
  max_history = min(nlalg.max_history, nlalg.max_iter, i)

  resize!(γs, max_history)
  resize!(Δz₊s, max_history)
  if max_history > max_history_old
    for i in (max_history_old + 1):max_history
      Δz₊s[i] = zero(z₊)
    end
  end

  if max_history != max_history_old
    nlcache.Q = typeof(nlcache.Q)(undef, i, max_history)
    nlcache.R = typeof(nlcache.R)(undef, max_history, max_history)
  end

  nothing
end

function Base.resize!(nlcache::NLAndersonConstantCache, nlalg::NLAnderson, i::Int)
  @unpack γs, Δz₊s = nlcache

  # determine new maximum history
  max_history_old = length(Δz₊s)
  max_history = min(nlalg.max_history, nlalg.max_iter, i)

  resize!(γs, max_history)
  resize!(Δz₊s, max_history)

  if max_history != max_history_old
    nlcache.Q = typeof(nlcache.Q)(undef, i, max_history)
    nlcache.R = typeof(nlcache.R)(undef, max_history, max_history)
  end

  nothing
end
