const DISCRETE_INPLACE_DEFAULT = convert(DiscreteFunction{true},(du,u,p,t) -> du.=u)
const DISCRETE_OUTOFPLACE_DEFAULT = convert(DiscreteFunction{false},(u,p,t) -> u)

"""
$(TYPEDEF)

Defines a discrete problem.

# Fields

$(FIELDS)
"""
struct DiscreteProblem{uType,tType,isinplace,P,F,C} <: AbstractDiscreteProblem{uType,tType,isinplace}
  """The function in the map."""
  f::F
  """The initial condition."""
  u0::uType
  """The timespan for the problem."""
  tspan::tType
  """The parameter values of the function."""
  p::P
  """ A callback to be applied to every solver which uses the problem."""
  callback::C
  @add_kwonly function DiscreteProblem{iip}(f::AbstractDiscreteFunction{iip},
                                            u0,tspan::Tuple,p=nothing;
                                            callback = nothing) where {iip}
    _tspan = promote_tspan(tspan)
    new{typeof(u0),typeof(_tspan),isinplace(f,4),
        typeof(p),
        typeof(f),typeof(callback)}(f,u0,_tspan,p,callback)
  end

  function DiscreteProblem{iip}(u0::Nothing,tspan::Nothing,p=nothing;
                                callback = nothing) where {iip}
    if iip
      f = DISCRETE_INPLACE_DEFAULT
    else
      f = DISCRETE_OUTOFPLACE_DEFAULT
    end
    new{Nothing,Nothing,iip,typeof(p),
        typeof(f),typeof(callback)}(f,nothing,nothing,p,callback)
  end

  function DiscreteProblem{iip}(f,u0,tspan,p=nothing;kwargs...) where {iip}
    DiscreteProblem(convert(DiscreteFunction{iip},f),u0,tspan,p;kwargs...)
  end
end

"""
    DiscreteProblem{isinplace}(f,u0,tspan,p=nothing,callback=nothing)

Defines a discrete problem with the specified functions.
"""
function DiscreteProblem(f::AbstractDiscreteFunction,u0,tspan::Tuple,p=nothing;kwargs...)
  DiscreteProblem{isinplace(f)}(f,u0,tspan,p;kwargs...)
end

function DiscreteProblem(f,u0,tspan::Tuple,p=nothing;kwargs...)
  iip = isinplace(f,4)
  DiscreteProblem(convert(DiscreteFunction{iip},f),u0,tspan,p;kwargs...)
end

"""
$(SIGNATURES)

Define a discrete problem with the identity map.
"""
function DiscreteProblem(u0,tspan::Tuple,p::Tuple;kwargs...)
  iip = typeof(u0) <: AbstractArray
  if iip
    f = DISCRETE_INPLACE_DEFAULT
  else
    f = DISCRETE_OUTOFPLACE_DEFAULT
  end
  DiscreteProblem(f,u0,tspan,p;kwargs...)
end

function DiscreteProblem(u0,tspan::Tuple,p=nothing;kwargs...)
  iip = typeof(u0) <: AbstractArray
  if iip
    f = DISCRETE_INPLACE_DEFAULT
  else
    f = DISCRETE_OUTOFPLACE_DEFAULT
  end
  DiscreteProblem(f,u0,tspan,p;kwargs...)
end
