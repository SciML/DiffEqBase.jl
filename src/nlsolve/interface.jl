## accessors
function nlsolve_f end

get_status(nlsolver::AbstractNLSolver) = nlsolver.status

nlsolvefail(nlsolver::AbstractNLSolver) = nlsolvefail(get_status(nlsolver))
nlsolvefail(nlstatus::NLStatus) = Int8(nlstatus) < 0

get_cache(nlsolver::AbstractNLSolver) = nlsolver.cache

get_new_W(nlsolver::AbstractNLSolver)::Bool = get_new_W(get_cache(nlsolver))
get_new_W(nlcache::AbstractNLSolverCache)::Bool = nlcache.new_W

set_new_W!(nlsolver::AbstractNLSolver, val::Bool)::Bool = set_new_W!(get_cache(nlsolver), val)
set_new_W!(nlcache::AbstractNLSolverCache, val::Bool)::Bool = (nlcache.new_W = val; val)

get_W(nlsolver::AbstractNLSolver) = get_W(get_cache(nlsolver))
get_W(nlcache::AbstractNLSolverCache) = nlcache.W

set_W!(nlsolver::AbstractNLSolver, W) = set_W!(get_cache(nlsolver), W)
set_W!(nlcache::AbstractNLSolverCache, W) = (nlcache.W = W; nothing)

get_W_dt(nlsolver::AbstractNLSolver) = get_W_dt(get_cache(nlsolver))
get_W_dt(nlcache::AbstractNLSolverCache) = nlcache.W_dt

set_W_dt!(nlsolver::AbstractNLSolver, W_dt) = set_W_dt!(get_cache(nlsolver), W_dt)
set_W_dt!(nlcache::AbstractNLSolverCache, W_dt) = (nlcache.W_dt = W_dt; nothing)

get_linsolve(nlsolver::AbstractNLSolver) = get_linsolve(get_cache(nlsolver))
get_linsolve(nlcache::AbstractNLSolverCache) = nlcache.linsolve

du_cache(nlsolver::AbstractNLSolver) = du_cache(get_cache(nlsolver))
du_cache(::AbstractNLSolverCache) = nothing

function get_nlsolver(integrator::DEIntegrator)
  isdefined(integrator.cache, :nlsolver) || return
  
  integrator.cache.nlsolver
end

## traits

isnewton(::AbstractNLSolver) = false

## build

function build_uf end
function build_jac_config end
function build_J_W end

## resize

function resize_jac_config! end
function resize_J! end
function resize_W! end

function resize_nlsolver!(integrator::DEIntegrator, i::Int)
  nlsolver = get_nlsolver(integrator)
 
  nlsolver === nothing && return

  if nlsolver isa AbstractArray
    for idx in eachindex(nlsolver)
      resize!(nlsolver[idx], integrator, i)
    end
  else
    resize!(nlsolver, integrator, i)
  end

  nothing
end

function Base.resize!(nlsolver::AbstractNLSolver, integrator, i::Int)
  @unpack z,zprev,tmp = nlsolver

  resize!(z, i)
  resize!(zprev, i)
  resize!(tmp, i)

  resize!(get_cache(nlsolver), nlsolver, integrator, i)
end

## default: dispatch only on the cache
Base.resize!(cache::AbstractNLSolverCache, nlsolver, integrator, i::Int) =
  Base.resize!(cache, i)