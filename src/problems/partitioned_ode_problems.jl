type PartitionedODEProblem{uType,tType,isinplace,F,C} <: AbstractPartitionedODEProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
end

function PartitionedODEProblem(f,u0,tspan; iip = isinplace(f[1],3),callback=nothing)
  @assert typeof(f) <: Tuple
  @assert typeof(u0) <: Tuple
  @assert length(f) == length(u0)
  PartitionedODEProblem{typeof(u0),eltype(tspan),iip,typeof(f),typeof(callback)}(f,u0,tspan,callback)
end

type SecondOrderODEProblem{uType,tType,isinplace,F,C} <: AbstractSecondOrderODEProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
end

function SecondOrderODEProblem(f,u0,tspan; iip = isinplace(f,4),callback=nothing)
  @assert typeof(u0) <: Tuple
  @assert length(u0) == 2
  SecondOrderODEProblem{typeof(u0),eltype(tspan),iip,typeof(f),typeof(callback)}(f,u0,tspan,callback)
end
