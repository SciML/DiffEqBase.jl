get_status(nlsolver::NLSolver) = nlsolver.status

nlsolvefail(nlsolver) = nlsolvefail(get_status(nlsolver))
nlsolvefail(nlstatus::NLStatus) = Int8(nlstatus) < 0

isnewton(::NLSolver) = false
isnewton(::NLSolver{<:NLNewton}) = true

get_new_W(nlsolver::NLSolver)::Bool = get_new_W(nlsolver.cache)
get_new_W(nlcache::Union{NLNewtonCache,NLNewtonConstantCache})::Bool = nlcache.new_W

set_new_W!(nlsolver::NLSolver, val::Bool)::Bool = set_new_W!(nlsolver.cache, val)
set_new_W!(nlcache::Union{NLNewtonCache,NLNewtonConstantCache}, val::Bool)::Bool = (nlcache.new_W = val; val)

get_W(nlsolver::NLSolver) = get_W(nlsolver.cache)
get_W(nlcache::Union{NLNewtonCache,NLNewtonConstantCache}) = nlcache.W

set_W!(nlsolver::NLSolver, W) = set_W!(nlsolver.cache, W)
set_W!(nlcache::Union{NLNewtonCache,NLNewtonConstantCache}, W) = (nlcache.W = W; W)

get_W_dt(nlsolver::NLSolver) = get_W_dt(nlsolver.cache)
get_W_dt(nlcache::Union{NLNewtonCache,NLNewtonConstantCache}) = nlcache.W_dt

set_W_dt!(nlsolver::NLSolver, W_dt) = set_W_dt!(nlsolver.cache, W_dt)
set_W_dt!(nlcache::Union{NLNewtonCache,NLNewtonConstantCache}, W_dt) = (nlcache.W_dt = W_dt; nothing)

get_linsolve(nlsolver::NLSolver) = get_linsolve(nlsolver.cache)
get_linsolve(nlcache::NLNewtonCache) = nlcache.linsolve
