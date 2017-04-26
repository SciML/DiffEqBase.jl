# Mu' = f
type ODEProblem{uType,tType,isinplace,F,C,MM} <: AbstractODEProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
  mass_matrix::MM
end

function ODEProblem(f,u0,tspan; iip = isinplace(f,3),callback=nothing,mass_matrix=I)
  ODEProblem{typeof(u0),promote_type(map(typeof,tspan)...),iip,typeof(f),typeof(callback),typeof(mass_matrix)}(f,u0,tspan,callback,mass_matrix)
end

# Mu' = f[1] + f[2] + ...
type SplitODEProblem{uType,tType,isinplace,F,C,MM} <: AbstractODEProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
  mass_matrix::MM
end

function SplitODEProblem(f,u0,tspan; iip = isinplace(f[2],3),callback=nothing,mass_matrix=I)
  @assert typeof(f) <: Tuple
  SplitODEProblem{typeof(u0),promote_type(map(typeof,tspan)...),iip,typeof(f),typeof(callback),typeof(mass_matrix)}(f,u0,tspan,callback,mass_matrix)
end

# M[i]*u[i]' = f[i]
type PartitionedODEProblem{uType,tType,isinplace,F,C,MM} <: AbstractODEProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
  mass_matrix::MM
end

function PartitionedODEProblem(f,u0,tspan; iip = isinplace(f[1],3),callback=nothing,mass_matrix=nothing)
  @assert typeof(f) <: Tuple
  @assert typeof(u0) <: Tuple
  if mass_matrix == nothing
    _mm = ((I for i in 1:length(f))...)
  end
  PartitionedODEProblem{typeof(u0),promote_type(map(typeof,tspan)...),iip,typeof(f),typeof(callback),typeof(_mm)}(f,u0,tspan,callback,_mm)
end

# u'' = f(t,u,du,ddu)
function SecondOrderODEProblem(f,u0,du0,tspan; iip = isinplace(f,4),callback=nothing,mass_matrix=I)
  if iip
    f1 = function (t,u,v,du)
      du .= v
    end
  else
    f1 = function (t,u,v)
      v
    end
  end
  _f = (f1,f)
  _u0 = (u0,du0)
  _mass_matrix = (mass_matrix,I)
  PartitionedODEProblem{typeof(_u0),promote_type(map(typeof,tspan)...),iip,typeof(_f),typeof(callback),typeof(_mass_matrix)}(_f,_u0,tspan,callback,_mass_matrix)
end
