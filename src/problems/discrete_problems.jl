const DISCRETE_INPLACE_DEFAULT = (t,u,du) -> fill!(du,zero(eltype(u)))
const DISCRETE_OUTOFPLACE_DEFAULT = (t,u) -> u

type DiscreteProblem{uType,tType,isinplace,F,C} <: AbstractDiscreteProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
  function DiscreteProblem{isinplace}(f,u0,tspan;callback = nothing) where {isinplace}
    new{typeof(u0),promote_type(map(typeof,tspan)...),isinplace,
        typeof(f),typeof(callback)}(f,u0,tspan,callback)
  end
end

function DiscreteProblem(f,u0,tspan;kwargs...)
  iip = isinplace(f,3)
  DiscreteProblem{iip}(f,u0,tspan;kwargs...)
end

function DiscreteProblem(u0,tspan::Tuple;kwargs...)
  iip = typeof(u0) <: AbstractArray
  if iip
    f = DISCRETE_INPLACE_DEFAULT
  else
    f = DISCRETE_OUTOFPLACE_DEFAULT
  end
  DiscreteProblem{iip}(f,u0,tspan;kwargs...)
end
