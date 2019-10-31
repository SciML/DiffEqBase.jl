# Mu' = f
"""
$(TYPEDEF)

Defines a steady state problem.

# Fields

$(FIELDS)
"""
struct SteadyStateProblem{uType,isinplace,P,F,K} <: AbstractSteadyStateProblem{uType,isinplace}
  """f: The function in the ODE."""
  f::F
  """The initial guess for the steady state."""
  u0::uType
  """Parameter values for the ODE function."""
  p::P
  kwargs::K
  @add_kwonly function SteadyStateProblem{iip}(f::AbstractODEFunction{iip},
                                               u0,p=NullParameters();
                                               kwargs...) where {iip}
    new{typeof(u0),isinplace(f),typeof(p),typeof(f),typeof(kwargs)}(f,u0,p,kwargs)
  end

  """
  $(SIGNATURES)

  Define a steady state problem using the given function.
  `isinplace` optionally sets whether the function is inplace or not.
  This is determined automatically, but not inferred.
  """
  function SteadyStateProblem{iip}(f,u0,p=NullParameters()) where iip
    SteadyStateProblem(ODEFunction{iip}(f),u0,p)
  end
end

"""
$(SIGNATURES)

Define a steady state problem using an instance of
[`AbstractODEFunction`](@ref DiffEqBase.AbstractODEFunction).
"""
function SteadyStateProblem(f::AbstractODEFunction,u0,p=NullParameters();kwargs...)
  SteadyStateProblem{isinplace(f)}(f,u0,p;kwargs...)
end

function SteadyStateProblem(f,u0,p=NullParameters();kwargs...)
  SteadyStateProblem(ODEFunction(f),u0,p;kwargs...)
end

"""
$(SIGNATURES)

Define a steady state problem from a standard ODE problem.
"""
SteadyStateProblem(prob::AbstractODEProblem) =
      SteadyStateProblem{isinplace(prob)}(prob.f,prob.u0,prob.p)
