# Mu' = f
type ODEProblem{uType,tType,isinplace,F,C,MM} <: AbstractODEProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
  mm::MM
end

function ODEProblem(f,u0,tspan; iip = isinplace(f,3),callback=nothing,mm=I)
  ODEProblem{typeof(u0),eltype(tspan),iip,typeof(f),typeof(callback),typeof(mm)}(f,u0,tspan,callback,mm)
end

# Mu' = f[1] + f[2] + ...
type SplitODEProblem{uType,tType,isinplace,F,C,MM} <: AbstractODEProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
  mm::MM
end

function SplitODEProblem(f,u0,tspan; iip = isinplace(f[2],3),callback=nothing,mm=I)
  @assert typeof(f) <: Tuple
  SplitODEProblem{typeof(u0),eltype(tspan),iip,typeof(f),typeof(callback),typeof(mm)}(f,u0,tspan,callback,mm)
end

# M[i]*u[i]' = f[i]
type PartitionedODEProblem{uType,tType,isinplace,F,C,MM} <: AbstractODEProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
  mm::MM
end

function PartitionedODEProblem(f,u0,tspan; iip = isinplace(f[1],3),callback=nothing,mm=nothing)
  @assert typeof(f) <: Tuple
  @assert typeof(u0) <: Tuple
  @assert length(f) == length(u0)
  if mm == nothing
    _mm = ((I for i in 1:length(f))...)
  end
  PartitionedODEProblem{typeof(u0),eltype(tspan),iip,typeof(f),typeof(callback),typeof(_mm)}(f,u0,tspan,callback,_mm)
end

# u'' = f(t,u,du,ddu)
type SecondOrderODEProblem{uType,tType,isinplace,F,C,MM} <: AbstractSecondOrderODEProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
  mm::MM
end

function SecondOrderODEProblem(f,u0,tspan; iip = isinplace(f,4),callback=nothing,mm=I)
  @assert typeof(u0) <: Tuple
  @assert length(u0) == 2
  SecondOrderODEProblem{typeof(u0),eltype(tspan),iip,typeof(f),typeof(callback),typeof(mm)}(f,u0,tspan,callback,mm)
end
