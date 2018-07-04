const DISCRETE_INPLACE_DEFAULT = (du,u,p,t) -> du.=u
const DISCRETE_OUTOFPLACE_DEFAULT = (u,p,t) -> u

struct DiscreteProblem{uType,tType,isinplace,P,F,C} <: AbstractDiscreteProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  p::P
  callback::C
  @add_kwonly function DiscreteProblem{iip}(f,u0,tspan,p=nothing;
           callback = nothing) where {iip}
    _tspan = promote_tspan(tspan)
    new{typeof(u0),typeof(_tspan),iip,
        typeof(p),
        typeof(f),typeof(callback)}(f,u0,_tspan,p,callback)
  end
end

function DiscreteProblem(f,u0,tspan::Tuple,p::Tuple;kwargs...)
  iip = isinplace(f,4)
  DiscreteProblem{iip}(f,u0,tspan,p;kwargs...)
end

function DiscreteProblem(f,u0,tspan::Tuple,p=nothing;kwargs...)
  iip = isinplace(f,4)
  DiscreteProblem{iip}(f,u0,tspan,p;kwargs...)
end

function DiscreteProblem(u0,tspan::Tuple,p::Tuple;kwargs...)
    iip = typeof(u0) <: AbstractArray
    if iip
      f = DISCRETE_INPLACE_DEFAULT
    else
      f = DISCRETE_OUTOFPLACE_DEFAULT
    end
    DiscreteProblem{iip}(f,u0,tspan,p;kwargs...)
end

function DiscreteProblem(u0,tspan::Tuple,p=nothing;kwargs...)
  iip = typeof(u0) <: AbstractArray
  if iip
    f = DISCRETE_INPLACE_DEFAULT
  else
    f = DISCRETE_OUTOFPLACE_DEFAULT
  end
  DiscreteProblem{iip}(f,u0,tspan,p;kwargs...)
end
