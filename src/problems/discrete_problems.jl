const DISCRETE_INPLACE_DEFAULT = (t,u,du) -> fill!(du,zero(eltype(u)))
const DISCRETE_OUTOFPLACE_DEFAULT = (t,u) -> u

type DiscreteProblem{uType,tType,isinplace,F,C} <: AbstractDiscreteProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
end

function DiscreteProblem(f,u0,tspan; iip = isinplace(f,3),callback = nothing)
  DiscreteProblem{typeof(u0),eltype(tspan),iip,typeof(f),typeof(callback)}(f,u0,tspan,callback)
end

function DiscreteProblem(u0,tspan::Tuple; callback = nothing)
  iip = typeof(u0) <: AbstractArray
  if iip
    f = DISCRETE_INPLACE_DEFAULT
  else
    f = DISCRETE_OUTOFPLACE_DEFAULT
  end
  DiscreteProblem{typeof(u0),eltype(tspan),iip,typeof(f),typeof(callback)}(f,u0,tspan,callback)
end
