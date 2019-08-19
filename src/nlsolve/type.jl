## algorithms

abstract type AbstractNLSolverAlgorithm end

struct NLFunctional{K,C} <: AbstractNLSolverAlgorithm
  κ::K
  fast_convergence_cutoff::C
  maxiters::Int
end

NLFunctional(; κ=1//100, maxiters=10, fast_convergence_cutoff=1//5) = NLFunctional(κ, fast_convergence_cutoff, maxiters)

struct NLAnderson{K,D,C} <: AbstractNLSolverAlgorithm
  κ::K
  fast_convergence_cutoff::C
  maxiters::Int
  max_history::Int
  aa_start::Int
  droptol::D
end

NLAnderson(; κ=1//100, maxiters=10, max_history::Int=10, aa_start::Int=1, droptol=1e10, fast_convergence_cutoff=1//5) =
  NLAnderson(κ, fast_convergence_cutoff, maxiters, max_history, aa_start, droptol)

struct NLNewton{K,C1,C2} <: AbstractNLSolverAlgorithm
  κ::K
  maxiters::Int
  fast_convergence_cutoff::C1
  new_W_dt_cutoff::C2
end

NLNewton(; κ=1//100, maxiters=10, fast_convergence_cutoff=1//5, new_W_dt_cutoff=1//5) =
  NLNewton(κ, maxiters, fast_convergence_cutoff, new_W_dt_cutoff)

## status

@enum NLStatus::Int8 begin
  FastConvergence     = 2
  Convergence         = 1
  SlowConvergence     = 0
  VerySlowConvergence = -1
  Divergence          = -2
  MaxItersReached      = -3
end

## NLsolver

abstract type AbstractNLSolver end

mutable struct NLSolver{algType<:AbstractNLSolverAlgorithm,IIP,uType,uTolType,tTypeNoUnits,C} <: AbstractNLSolver
  z::uType
  zprev::uType
  tmp::uType
  γ::uTolType
  c::tTypeNoUnits
  alg::algType
  κ::uTolType
  η::uTolType
  fast_convergence_cutoff::uTolType
  iter::Int
  maxiters::Int
  status::NLStatus
  cache::C
end

## caches

abstract type AbstractNLSolverCache end

mutable struct NLFunctionalCache{uType,tType,rateType,uNoUnitsType} <: AbstractNLSolverCache
  dz::uType
  tstep::tType
  k::rateType
  atmp::uNoUnitsType
end

mutable struct NLFunctionalConstantCache{uType,tType} <: AbstractNLSolverCache
  dz::uType
  tstep::tType
end

mutable struct NLAndersonCache{uType,tType,rateType,uNoUnitsType,uEltypeNoUnits,D} <: AbstractNLSolverCache
  dz::uType
  tstep::tType
  k::rateType
  atmp::uNoUnitsType
  gprev::uType
  dzprev::uType
  Δgs::Vector{uType}
  Q::Matrix{uEltypeNoUnits}
  R::Matrix{uEltypeNoUnits}
  γs::Vector{uEltypeNoUnits}
  history::Int
  droptol::D
end

mutable struct NLAndersonConstantCache{uType,tType,uEltypeNoUnits,D} <: AbstractNLSolverCache
  dz::uType
  tstep::tType
  gprev::uType
  dzprev::uType
  Δgs::Vector{uType}
  Q::Matrix{uEltypeNoUnits}
  R::Matrix{uEltypeNoUnits}
  γs::Vector{uEltypeNoUnits}
  history::Int
  droptol::D
end

mutable struct NLNewtonCache{uType,tType,rateType,uNoUnitsType,J,W,du1Type,ufType,jcType,lsType,G} <: AbstractNLSolverCache
  dz::uType
  tstep::tType
  k::rateType
  atmp::uNoUnitsType
  J::J
  W::W
  new_W::Bool
  W_dt::tType
  du1::du1Type
  uf::ufType
  jac_config::jcType
  linsolve::lsType
  weight::uType
  invγdt::G
  new_W_dt_cutoff::tType
end

mutable struct NLNewtonConstantCache{uType,tType,J,W,ufType,G} <: AbstractNLSolverCache
  dz::uType
  tstep::tType
  J::J
  W::W
  new_W::Bool
  W_dt::tType
  uf::ufType
  invγdt::G
  new_W_dt_cutoff::tType
end