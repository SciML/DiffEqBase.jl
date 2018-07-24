# Mu' = f
struct SteadyStateProblem{uType,isinplace,P,F} <: AbstractSteadyStateProblem{uType,isinplace}
  f::F
  u0::uType
  p::P
  @add_kwonly function SteadyStateProblem{iip}(f,u0,p=nothing) where iip
    new{typeof(u0),iip,typeof(p),
        typeof(f)}(f,u0,p)
  end
end

function SteadyStateProblem(f,u0,p=nothing;kwargs...)
  iip = isinplace(f,3)
  SteadyStateProblem{iip}(f,u0,p;kwargs...)
end

SteadyStateProblem(prob::AbstractODEProblem) =
      SteadyStateProblem{isinplace(prob)}(prob.f,prob.u0)
