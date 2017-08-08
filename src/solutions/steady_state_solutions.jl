mutable struct SteadyStateSolution{T,N,uType,R,P,A} <: AbstractSteadyStateSolution{T,N}
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
  if typeof(prob.u0) <: Tuple
    N = length((size(ArrayPartition(prob.u0))...))
  else
    N = length((size(prob.u0)...))
  end

  SteadyStateSolution{T,N,typeof(u),typeof(resid),typeof(prob),typeof(alg)}(u,resid,prob,alg,retcode)
end
