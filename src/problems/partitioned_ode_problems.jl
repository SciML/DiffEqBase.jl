type PartitionedODEProblem{uType,tType,isinplace,F,C} <: AbstractPartitionedODEProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
end

function PartitionedODEProblem(f,u0,tspan; iip = isinplace(f[1],3),callback=CallbackSet())
  @assert typeof(f) <: Tuple
  @assert typeof(u0) <: Tuple
  @assert length(f) == length(u0)
  PartitionedODEProblem{typeof(u0),eltype(tspan),iip,typeof(f),typeof(callback)}(f,u0,tspan,callback)
end
