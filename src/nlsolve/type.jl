abstract type AbstractNLSolverAlgorithm end
abstract type AbstractNLSolverCache end

@enum NLStatus::Int8 begin
  FastConvergence     = 2
  Convergence         = 1
  SlowConvergence     = 0
  VerySlowConvergence = -1
  Divergence          = -2
end

# solver

mutable struct NLSolver{iip,uType,rateType,uTolType,kType,gType,cType,du1Type,ufType,jcType,lsType,C1,C<:AbstractNLSolverCache}
  z::uType
  dz::uType
  tmp::uType
  ztmp::uType
  k::rateType
  ηold::uTolType
  κ::kType
  γ::gType
  c::cType
  max_iter::Int
  nl_iters::Int
  status::NLStatus
  fast_convergence_cutoff::C1
  du1::du1Type
  uf::ufType
  jac_config::jcType
  linsolve::lsType
  weight::uType
  cache::C
end

nlsolve_default(_, ::Val{:κ}) = 1//100 # takes `integrator.alg`
nlsolve_default(_, ::Val{:max_iter}) = 10
nlsolve_default(_, ::Val{:fast_convergence_cutoff}) = 1//5
nlsolve_default(_, ::Val{:new_W_dt_cutoff}) = 1//5
nlsolve_default(_, ::Val{:max_history}) = 10
nlsolve_default(_, ::Val{:aa_start}) = 1
nlsolve_default(_, ::Val{:droptol}) = Some(nothing)

# algorithms

struct NLFunctional{K,C,M} <: AbstractNLSolverAlgorithm
  κ::K
  fast_convergence_cutoff::C
  max_iter::M
end

NLFunctional(; κ=nothing, max_iter=nothing, fast_convergence_cutoff=nothing) = NLFunctional(κ, fast_convergence_cutoff, max_iter)

function handle_defaults(alg, nlalg::NLFunctional)
  @set! nlalg.κ = something(nlalg.κ, nlsolve_default(alg, Val(:κ)))
  @set! nlalg.max_iter = something(nlalg.max_iter, nlsolve_default(alg, Val(:max_iter)))
  @set! nlalg.fast_convergence_cutoff = something(nlalg.fast_convergence_cutoff, nlsolve_default(alg, Val(:fast_convergence_cutoff)))
  return nlalg
end

struct NLAnderson{K,D,C,MI,MH,AS} <: AbstractNLSolverAlgorithm
  κ::K
  fast_convergence_cutoff::C
  max_iter::MI
  max_history::MH
  aa_start::AS
  droptol::D
end

NLAnderson(; κ=nothing, max_iter=nothing, max_history=nothing, aa_start=nothing, droptol=nothing, fast_convergence_cutoff=nothing) =
  NLAnderson(κ, fast_convergence_cutoff, max_iter, max_history, aa_start, droptol)

function handle_defaults(alg, nlalg::NLAnderson)
  @set! nlalg.κ = something(nlalg.κ, nlsolve_default(alg, Val(:κ)))
  @set! nlalg.max_iter = something(nlalg.max_iter, nlsolve_default(alg, Val(:max_iter)))
  @set! nlalg.fast_convergence_cutoff = something(nlalg.fast_convergence_cutoff, nlsolve_default(alg, Val(:fast_convergence_cutoff)))
  @set! nlalg.max_history = something(nlalg.max_history, nlsolve_default(alg, Val(:max_history)))
  @set! nlalg.aa_start = something(nlalg.aa_start, nlsolve_default(alg, Val(:aa_start)))
  @set! nlalg.droptol = something(nlalg.droptol, nlsolve_default(alg, Val(:droptol)))
  return nlalg
end

struct NLNewton{K,M,C1,C2} <: AbstractNLSolverAlgorithm
  κ::K
  max_iter::M
  fast_convergence_cutoff::C1
  new_W_dt_cutoff::C2
end

NLNewton(; κ=nothing, max_iter=nothing, fast_convergence_cutoff=nothing, new_W_dt_cutoff=nothing) = NLNewton(κ, max_iter, fast_convergence_cutoff, new_W_dt_cutoff)

function handle_defaults(alg, nlalg::NLNewton)
  @set! nlalg.κ = something(nlalg.κ, nlsolve_default(alg, Val(:κ)))
  @set! nlalg.max_iter = something(nlalg.max_iter, nlsolve_default(alg, Val(:max_iter)))
  @set! nlalg.fast_convergence_cutoff = something(nlalg.fast_convergence_cutoff, nlsolve_default(alg, Val(:fast_convergence_cutoff)))
  @set! nlalg.new_W_dt_cutoff = something(nlalg.new_W_dt_cutoff, nlsolve_default(alg, Val(:new_W_dt_cutoff)))
  return nlalg
end

# caches

mutable struct NLNewtonCache{W,J,T,C} <: AbstractNLSolverCache
  new_W::Bool
  W::W
  J::J
  W_dt::T
  new_W_dt_cutoff::C
end

mutable struct NLNewtonConstantCache{W,J,C} <: AbstractNLSolverCache
  W::W
  J::J
  new_W_dt_cutoff::C
end

struct NLFunctionalCache{uType} <: AbstractNLSolverCache
  z₊::uType
end

struct NLFunctionalConstantCache <: AbstractNLSolverCache end

mutable struct NLAndersonCache{uType,gsType,QType,RType,gType,D} <: AbstractNLSolverCache
  z₊::uType
  dzold::uType
  z₊old::uType
  Δz₊s::gsType
  Q::QType
  R::RType
  γs::gType
  aa_start::Int
  droptol::D
end

mutable struct NLAndersonConstantCache{gsType,QType,RType,gType,D} <: AbstractNLSolverCache
  Δz₊s::gsType
  Q::QType
  R::RType
  γs::gType
  aa_start::Int
  droptol::D
end
