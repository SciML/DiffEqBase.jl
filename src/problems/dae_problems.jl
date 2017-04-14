# f(t,u,du,res) = 0
type DAEProblem{uType,duType,tType,isinplace,F,C} <: AbstractDAEProblem{uType,duType,tType,isinplace}
  f::F
  u0::uType
  du0::duType
  tspan::Tuple{tType,tType}
  callback::C
end

function DAEProblem(f,u0,du0,tspan; iip = isinplace(f,4), callback = nothing,mm=I,diff_vars=nothing)
  DAEProblem{typeof(u0),typeof(du0),eltype(tspan),iip,typeof(f),typeof(callback)}(f,u0,du0,tspan,callback)
end

# 0 = f[1](t,u,du,res) + f[2](t,u,du,res) + ...
type SplitDAEProblem{uType,duType,tType,isinplace,F,C} <: AbstractDAEProblem{uType,duType,tType,isinplace}
  f::F
  u0::uType
  du0::duType
  tspan::Tuple{tType,tType}
  callback::C
end

function SplitDAEProblem(f,u0,du0,tspan; iip = isinplace(f,4), callback = nothing,mm=I)
  @assert typeof(f) <: Tuple
  SplitDAEProblem{typeof(u0),typeof(du0),eltype(tspan),iip,typeof(f),typeof(callback)}(f,u0,du0,tspan,callback)
end

# Implicit problem for each variable in the partition
# 0 = f[i](t,u...,du[1],res), 0 = f[2](t,u...,du[2],res), ...
type PartitionedDAEProblem{uType,duType,tType,isinplace,F,C} <: AbstractDAEProblem{uType,duType,tType,isinplace}
  f::F
  u0::uType
  du0::duType
  tspan::Tuple{tType,tType}
  callback::C
end

function PartitionedDAEProblem(f,u0,du0,tspan; iip = isinplace(f,4), callback = nothing,mm=I)
  @assert typeof(f) <: Tuple
  SplitDAEProblem{typeof(u0),typeof(du0),eltype(tspan),iip,typeof(f),typeof(callback)}(f,u0,du0,tspan,callback)
end
