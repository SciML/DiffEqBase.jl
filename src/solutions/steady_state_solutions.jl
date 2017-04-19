type SteadyStateSolution{uType,P,A} <: AbstractSteadyStateSolution
  u::uType
  prob::P
  alg::A
  retcode::Symbol
end

function build_solution(prob::AbstractSteadyStateProblem,
                        alg,u;calculate_error = true,
                        retcode = :Default, kwargs...)
  SteadyStateSolution(u,prob,alg,retcode)
end
