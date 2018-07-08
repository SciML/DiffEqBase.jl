const DISCRETE_INPLACE_DEFAULT = convert(DiscreteFunction{true},(du,u,p,t) -> du.=u)
const DISCRETE_OUTOFPLACE_DEFAULT = convert(DiscreteFunction{false},(u,p,t) -> u)

struct DiscreteProblem{uType,tType,isinplace,P,F,C} <: AbstractDiscreteProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::tType
  p::P
  callback::C
  @add_kwonly function DiscreteProblem(f::AbstractDiscreteFunction,
                                            u0,tspan::Tuple,p=nothing;
                                            callback = nothing)
    _tspan = promote_tspan(tspan)
    new{typeof(u0),typeof(_tspan),isinplace(f,4),
        typeof(p),
        typeof(f),typeof(callback)}(f,u0,_tspan,p,callback)
  end

  function DiscreteProblem{iip}(f,u0,tspan,p=nothing;kwargs...) where {iip}
    DiscreteProblem(convert(DiscreteFunction{iip},f),u0,tspan,p;kwargs...)
  end
end

function DiscreteProblem(f,u0,tspan::Tuple,p;kwargs...)
  iip = isinplace(f,4)
  DiscreteProblem(convert(DiscreteFunction{iip},f),u0,tspan,p;kwargs...)
end

function DiscreteProblem(f,u0,tspan::Tuple,p=nothing;kwargs...)
  iip = isinplace(f,4)
  DiscreteProblem(convert(DiscreteFunction{iip},f),u0,tspan,p;kwargs...)
end

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
