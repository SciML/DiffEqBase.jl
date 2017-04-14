# u'=f(t,u,v,du) with 0 = g(t,u,v)
type ConstrainedODEProblem{uType,duType,tType,isinplace,F,G,C,MM} <: AbstractDAEProblem{uType,duType,tType,isinplace}
  f::F
  g::G
  u0::uType
  v0::duType
  tspan::Tuple{tType,tType}
  callback::C
  mm::MM
end

function ConstrainedODEProblem(f,g,u0,v0,tspan; iip = isinplace(f,4), callback = nothing,mm=I)
  ConstrainedODEProblem{typeof(u0),typeof(v0),eltype(tspan),iip,typeof(f),typeof(g),typeof(callback),typeof(mm)}(f,g,u0,v0,tspan,callback,mm)
end

# u' = f[1](t,u,v,du) + f[2](t,u,v,du) + ... and a single 0 = g(t,u,v)
type SplitConstrainedODEProblem{uType,duType,tType,isinplace,F,G,C,MM} <: AbstractDAEProblem{uType,duType,tType,isinplace}
  f::F
  g::G
  u0::uType
  v0::duType
  tspan::Tuple{tType,tType}
  callback::C
  mm::MM
end

function SplitConstrainedODEProblem(f,g,u0,v0,tspan; iip = isinplace(f,4), callback = nothing,mm=I)
  @assert typeof(f) <: Tuple
  @assert typeof(g) <: Tuple
  SplitConstrainedODEProblem{typeof(u0),typeof(v0),eltype(tspan),iip,typeof(f),typeof(g),typeof(callback),typeof(_mm)}(f,g,u0,v0,tspan,callback,mm)
end

# Tuple of u[i]'=f[1] and a single 0 = g(t,u,v)
type PartitionedConstrainedODEProblem{uType,duType,tType,isinplace,F,G,C,MM} <: AbstractDAEProblem{uType,duType,tType,isinplace}
  f::F
  g::G
  u0::uType
  du0::duType
  tspan::Tuple{tType,tType}
  callback::C
  mm::MM
end

function PartitionedConstrainedODEProblem(f,g,u0,du0,tspan; iip = isinplace(f,4), callback = nothing,mm=nothing)
  @assert typeof(f) <: Tuple
  @assert typeof(g) <: Tuple
  if mm == nothing
    _mm = ((I for i in 1:length(f))...)
  end
  PartitionedConstrainedODEProblem{typeof(u0),typeof(du0),eltype(tspan),iip,typeof(f),typeof(g),typeof(callback),typeof(_mm)}(f,g,u0,du0,tspan,callback,_mm)
end
