# Mu' = f
struct SteadyStateProblem{uType,isinplace,P,F} <: AbstractSteadyStateProblem{uType,isinplace}
  f::F
  u0::uType
  p::P
  @add_kwonly function SteadyStateProblem(f::AbstractODEFunction,u0,p=nothing)
    new{typeof(u0),isinplace(f),typeof(p),typeof(f)}(f,u0,p)
  end

  function SteadyStateProblem{iip}(f,u0,p=nothing) where iip
    SteadyStateProblem(ODEFunction{iip}(f),u0,p)
  end
end

function SteadyStateProblem(f,u0,p=nothing;kwargs...)
  SteadyStateProblem(ODEFunction(f),u0,p;kwargs...)
end

SteadyStateProblem(prob::AbstractODEProblem) =
      SteadyStateProblem{isinplace(prob)}(prob.f,prob.u0)
