type SDEProblem{uType,tType,isinplace,NP,F,F2,C,MM,ND} <: AbstractSDEProblem{uType,tType,isinplace,ND}
  f::F
  g::F2
  u0::uType
  tspan::Tuple{tType,tType}
  noise::NP
  callback::C
  mass_matrix::MM
  noise_rate_prototype::ND
end

function SDEProblem(f,g,u0,tspan;iip = isinplace(f,3),
        noise_rate_prototype = nothing,
        noise= nothing,
        callback = nothing,mass_matrix=I)

  SDEProblem{typeof(u0),eltype(tspan),iip,typeof(noise),typeof(f),typeof(g),typeof(callback),typeof(mass_matrix),typeof(noise_rate_prototype)}(
  f,g,u0,tspan,noise,callback,mass_matrix,noise_rate_prototype)
end

# Mu' = f[1] + f[2] + ...
type SplitSDEProblem{uType,tType,isinplace,F,C,MM} <: AbstractSDEProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
  mass_matrix::MM
end

function SplitSDEProblem(f,u0,tspan; iip = isinplace(f[2],3),callback=nothing,mass_matrix=I)
  @assert typeof(f) <: Tuple
  SplitSDEProblem{typeof(u0),eltype(tspan),iip,typeof(f),typeof(callback),typeof(mass_matrix)}(f,u0,tspan,callback,mass_matrix)
end

# M[i]*u[i]' = f[i]
type PartitionedSDEProblem{uType,tType,isinplace,F,C,MM} <: AbstractSDEProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
  mass_matrix::MM
end

function PartitionedSDEProblem(f,u0,tspan; iip = isinplace(f[1],3),callback=nothing,mass_matrix=nothing)
  @assert typeof(f) <: Tuple
  @assert typeof(u0) <: Tuple
  if mass_matrix == nothing
    _mm = ((I for i in 1:length(f))...)
  end
  PartitionedSDEProblem{typeof(u0),eltype(tspan),iip,typeof(f),typeof(callback),typeof(_mm)}(f,u0,tspan,callback,_mm)
end

# u'' = f(t,u,du,ddu)
function SecondOrderSDEProblem(f,u0,du0,tspan; iip = isinplace(f,4),callback=nothing,mass_matrix=I)
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
  PartitionedSDEProblem{typeof(_u0),eltype(tspan),iip,typeof(_f),typeof(callback),typeof(_mass_matrix)}(_f,_u0,tspan,callback,_mass_matrix)
end
