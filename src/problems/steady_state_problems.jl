# Mu' = f
mutable struct SteadyStateProblem{uType,isinplace,F,MM} <: AbstractSteadyStateProblem{uType,isinplace}
  f::F
  u0::uType
  mass_matrix::MM
  function SteadyStateProblem{iip}(f,u0;mass_matrix=I) where iip
    SteadyStateProblem{typeof(u0),iip,typeof(f),typeof(mass_matrix)}(f,u0,mass_matrix)
  end
end

function SteadyStateProblem(f,u0;kwargs...)
  iip = isinplace(f,3)
  SteadyStateProblem{iip}(f,u0;kwargs...)
end

SteadyStateProblem(prob::AbstractODEProblem) =
      SteadyStateProblem(prob.f,prob.u0,iip=isinplace(prob),mass_matrix=prob.mass_matrix)
