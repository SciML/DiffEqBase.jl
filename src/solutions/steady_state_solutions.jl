type SteadyStateSolution{uType,R,P,A} <: AbstractSteadyStateSolution
  u::uType
  resid::R
  prob::P
  alg::A
  retcode::Symbol
end

function build_solution(prob::AbstractSteadyStateProblem,
                        alg,u,resid;calculate_error = true,
                        retcode = :Default, kwargs...)
  SteadyStateSolution(u,resid,prob,alg,retcode)
end
