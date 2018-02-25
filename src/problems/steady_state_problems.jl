# Mu' = f
struct SteadyStateProblem{uType,isinplace,P,J,F,MM} <: AbstractSteadyStateProblem{uType,isinplace}
  f::F
  u0::uType
  p::P
  jac_prototype::J
  mass_matrix::MM
  @add_kwonly function SteadyStateProblem{iip}(f,u0,p=nothing;
                                   jac_prototype=nothing,
                                   mass_matrix=I) where iip
    new{typeof(u0),iip,typeof(p),typeof(jac_prototype),
        typeof(f),typeof(mass_matrix)}(f,u0,p,jac_prototype,mass_matrix)
  end
end

function SteadyStateProblem(f,u0,p=nothing;kwargs...)
  iip = isinplace(f,3)
  SteadyStateProblem{iip}(f,u0,p;kwargs...)
end

SteadyStateProblem(prob::AbstractODEProblem) =
      SteadyStateProblem{isinplace(prob)}(prob.f,prob.u0,mass_matrix=prob.mass_matrix)
