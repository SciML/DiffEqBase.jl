struct LinearSolution{T,N,uType,R,P,A} <: AbstractLinearSolution{T,N}
  u::uType
  resid::R
  prob::P
  alg::A
  retcode::Symbol
end

function build_solution(prob::AbstractLinearProblem,
                        alg,u,resid;calculate_error = true,
                        retcode = :Default, kwargs...)

  T = eltype(eltype(u))
  N = length((size(u)...,))

  LinearSolution{T,N,typeof(u),typeof(resid),typeof(prob),typeof(alg)}(u,resid,prob,alg,retcode)
end

struct NonlinearSolution{T,N,uType,R,P,A} <: AbstractNonlinearSolution{T,N}
  u::uType
  resid::R
  prob::P
  alg::A
  retcode::Symbol
end

function build_solution(prob::AbstractNonlinearProblem,
                        alg,u,resid;calculate_error = true,
                        retcode = :Default, kwargs...)

  T = eltype(eltype(u))
  N = length((size(u)...,))

  NonlinearSolution{T,N,typeof(u),typeof(resid),typeof(prob),typeof(alg)}(u,resid,prob,alg,retcode)
end

struct QuadratureSolution{T,N,uType,R,P,A,C} <: AbstractQuadratureSolution{T,N}
  u::uType
  resid::R
  prob::P
  alg::A
  retcode::Symbol
  chi::C
end

function build_solution(prob::AbstractQuadratureProblem,
                        alg,u,resid;calculate_error = true,
                        chi = nothing,
                        retcode = :Default, kwargs...)

  T = eltype(eltype(u))
  N = length((size(u)...,))

  QuadratureSolution{T,N,typeof(u),typeof(resid),
                     typeof(prob),typeof(alg),typeof(chi)}(
                     u,resid,prob,alg,retcode,chi)
end
