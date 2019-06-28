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

get_status(nlsolver::NLSolver) = nlsolver.status
nlsolvefail(nlsolver::NLSolver) = nlsolvefail(get_status(nlsolver))
nlsolvefail(status::NLStatus) = Int8(status) < 0

isnewton(nlsolver::NLSolver) = isnewton(nlsolver.cache)
isnewton(nlcache::Union{NLNewtonCache,NLNewtonConstantCache}) = true
isnewton(nlcache::AbstractNLSolverCache) = false

set_new_W!(nlsolver::NLSolver, val::Bool)::Bool = set_new_W!(nlsolver.cache, val)
set_new_W!(nlcache::NLNewtonCache, val::Bool)::Bool = nlcache.new_W = val
set_new_W!(nlcache::AbstractNLSolverCache, val::Bool)::Bool = val

get_W(nlsolver::NLSolver) = get_W(nlsolver.cache)
get_W(nlcache::Union{NLNewtonCache,NLNewtonConstantCache}) = nlcache.W

set_W!(nlsolver::NLSolver, W) = set_W!(nlsolver.cache, W)
set_W!(nlcache::Union{NLNewtonCache,NLNewtonConstantCache}, W) = (nlcache.W = W; W)

set_W_dt!(nlsolver::NLSolver, W_dt) = set_W_dt!(nlsolver.cache, W_dt)
set_W_dt!(nlcache::NLNewtonCache, W_dt) = (nlcache.W_dt = W_dt; W_dt)
set_W_dt!(nlcache::NLNewtonConstantCache, W_dt) = W_dt

function nlsolve_f end

function iipnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @unpack κ, fast_convergence_cutoff = alg.nlsolve

  # define additional fields of cache of non-linear solver
  z = similar(u); dz = similar(u); tmp = similar(u); b = similar(u)
  k = zero(rate_prototype)

  uTolType = real(uBottomEltypeNoUnits)

  # create cache of non-linear solver
  if alg.nlsolve isa NLNewton
    nlcache = DiffEqBase.NLNewtonCache(true,W,dt,alg.nlsolve.new_W_dt_cutoff)
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

    nlcache = DiffEqBase.NLAndersonCache(z₊,dzold,z₊old,Δz₊s,Q,R,γs,alg.nlsolve.aa_start,alg.nlsolve.droptol)
  end

  # define additional fields of cache
  if alg.nlsolve isa NLNewton
    # TODO: check if the solver is iterative
    weight = similar(u)
  else
    weight = z
  end

  # create non-linear solver
  nlsolver = NLSolver{true,typeof(z),typeof(k),uTolType,typeof(κ),typeof(γ),typeof(c),typeof(fast_convergence_cutoff),typeof(nlcache)}(z,dz,tmp,b,k,one(uTolType),κ,γ,c,alg.nlsolve.max_iter,10000,Convergence,fast_convergence_cutoff,weight,nlcache)
end

DiffEqBase.@def getiipnlsolvefields begin
  @unpack z,dz,tmp,k = nlsolver
  b = nlsolver.ztmp
  fsalfirst = zero(rate_prototype)

  if alg.nlsolve isa NLNewton
    nf = nlsolve_f(f, alg)
    islin = f isa Union{ODEFunction,SplitFunction} && islinear(nf.f)
    if islin
      du1 = rate_prototype
      uf = nothing
      jac_config = nothing
      linsolve = alg.linsolve(Val{:init},nf,u)
    else
      du1 = zero(rate_prototype)
      # if the algorithm specializes on split problems the use `nf`
      uf = DiffEqDiffTools.UJacobianWrapper(nf,t,p)
      jac_config = build_jac_config(alg,nf,uf,du1,uprev,u,tmp,dz)
      linsolve = alg.linsolve(Val{:init},uf,u)
    end
  else
    du1 = rate_prototype
    uf = nothing
    jac_config = nothing
    linsolve = nothing
  end
end

function oopnlsolve(alg,u,uprev,p,t,dt,f,W,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,γ,c)
  @unpack κ, fast_convergence_cutoff = alg.nlsolve

  # define additional fields of cache of non-linear solver (all aliased)
  z = uprev; dz = z; tmp = z; b = z; k = rate_prototype

  uTolType = real(uBottomEltypeNoUnits)

  # create cache of non-linear solver
  if alg.nlsolve isa NLNewton
    nlcache = DiffEqBase.NLNewtonConstantCache(W,alg.nlsolve.new_W_dt_cutoff)
  elseif alg.nlsolve isa NLFunctional
    nlcache = DiffEqBase.NLFunctionalConstantCache()
  elseif alg.nlsolve isa NLAnderson
    max_history = min(alg.nlsolve.max_history, alg.nlsolve.max_iter, length(z))
    Δz₊s = Vector{typeof(z)}(undef, max_history)
    Q = Matrix{uEltypeNoUnits}(undef, length(z), max_history)
    R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
    γs = Vector{uEltypeNoUnits}(undef, max_history)

    nlcache = DiffEqBase.NLAndersonConstantCache(Δz₊s,Q,R,γs,alg.nlsolve.aa_start,alg.nlsolve.droptol)
  end

  # create non-linear solver
  nlsolver = NLSolver{false,typeof(z),typeof(k),uTolType,typeof(κ),typeof(γ),typeof(c),typeof(fast_convergence_cutoff),typeof(nlcache)}(z,dz,tmp,b,k,one(uTolType),κ,γ,c,alg.nlsolve.max_iter,10000,Convergence,fast_convergence_cutoff,z,nlcache)
end

DiffEqBase.@def getoopnlsolvefields begin
    if alg.nlsolve isa NLNewton
      nf = nlsolve_f(f, alg)
      # only use `nf` if the algorithm specializes on split eqs
      uf = DiffEqDiffTools.UDerivativeWrapper(nf,t,p)
    else
      uf = nothing
    end
end