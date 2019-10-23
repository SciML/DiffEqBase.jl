"""
$(TYPEDEF)

TODO
"""
struct SteadyStateSolution{T,N,uType,R,P,A} <: AbstractSteadyStateSolution{T,N}
  u::uType
  resid::R
  prob::P
  alg::A
  retcode::Symbol
end

function build_solution(prob::AbstractSteadyStateProblem,
                        alg,u,resid;calculate_error = true,
                        retcode = :Default, kwargs...)

  T = eltype(eltype(u))
  N = length((size(prob.u0)...,))

  SteadyStateSolution{T,N,typeof(u),typeof(resid),typeof(prob),typeof(alg)}(u,resid,prob,alg,retcode)
end
