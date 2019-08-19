## accessors
function nlsolve_f end

get_status(nlsolver::NLSolver) = nlsolver.status

nlsolvefail(nlsolver) = nlsolvefail(get_status(nlsolver))
nlsolvefail(nlstatus::NLStatus) = Int8(nlstatus) < 0

get_new_W(nlsolver::NLSolver)::Bool = get_new_W(nlsolver.cache)
set_new_W!(nlsolver::NLSolver, val::Bool)::Bool = set_new_W!(nlsolver.cache, val)

get_W(nlsolver::NLSolver) = get_W(nlsolver.cache)
set_W!(nlsolver::NLSolver, W) = set_W!(nlsolver.cache, W)

get_W_dt(nlsolver::NLSolver) = get_W_dt(nlsolver.cache)
set_W_dt!(nlsolver::NLSolver, W_dt) = set_W_dt!(nlsolver.cache, W_dt)

get_linsolve(nlsolver::NLSolver) = get_linsolve(nlsolver.cache)

## traits

isnewton(::NLSolver) = false

## build

function build_uf end
function build_jac_config end
function build_J_W end

## resize

function resize_jac_config! end
function resize_J! end
function resize_W! end

function resize_nlsolver!(integrator::DEIntegrator, i::Int)
  isdefined(integrator.cache, :nlsolver) || return

  nlsolver = integrator.cache.nlsolver
  if nlsolver isa AbstractArray
    for idx in eachindex(nlsolver)
      resize!(nlsolver[idx], integrator, i)
    end
  else
    resize!(nlsolver, integrator, i)
  end

  nothing
end

function Base.resize!(nlsolver::NLSolver, integrator, i::Int)
  @unpack z,zprev,tmp,cache = nlsolver

  resize!(z, i)
  resize!(zprev, i)
  resize!(tmp, i)

  resize!(cache, nlsolver, integrator, i)
end

## default: dispatch only on the cache
Base.resize!(cache::AbstractNLSolverCache, nlsolver, integrator, i::Int) =
  Base.resize!(cache, i)