const DISCRETE_INPLACE_DEFAULT = (t,u,du) -> fill!(du,zero(eltype(u)))
const DISCRETE_OUTOFPLACE_DEFAULT = (t,u) -> u

type DiscreteProblem{uType,tType,isinplace,F} <: AbstractDiscreteProblem{uType,tType,isinplace,F}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
end

type DiscreteTestProblem{uType,AType,tType,isinplace,F} <: AbstractDiscreteTestProblem{uType,tType,isinplace,F}
  f::F
  u0::uType
  analytic::AType
  tspan::Tuple{tType,tType}
end

function DiscreteProblem(f,u0,tspan; iip = isinplace(f,3))
  DiscreteProblem{typeof(u0),eltype(tspan),iip,typeof(f)}(f,u0,tspan)
end

function DiscreteTestProblem(f,u0,analytic,tspan::Tuple=(0.0,1.0); iip = isinplace(f,3))
  DiscreteTestProblem{typeof(u0),typeof(analytic),eltype(tspan),iip,typeof(f)}(f,u0,analytic,tspan)
end

function DiscreteProblem(u0,tspan::Tuple)
  iip = typeof(u0) <: AbstractArray
  if iip
    f = DISCRETE_INPLACE_DEFAULT
  else
    f = DISCRETE_OUTOFPLACE_DEFAULT
  end
  DiscreteProblem{typeof(u0),eltype(tspan),iip,typeof(f)}(f,u0,tspan)
end

function DiscreteTestProblem(u0,analytic,tspan::Tuple=(0.0,1.0))
  iip = typeof(u0) <: AbstractArray
  if iip
    f = DISCRETE_INPLACE_DEFAULT
  else
    f = DISCRETE_OUTOFPLACE_DEFAULT
  end
  DiscreteTestProblem{typeof(u0),typeof(analytic),eltype(tspan),iip,typeof(f)}(f,u0,analytic,tspan)
end
