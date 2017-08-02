# Mu' = f
mutable struct ODEProblem{uType,tType,isinplace,F,C,MM} <: AbstractODEProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
  mass_matrix::MM
  function ODEProblem{isinplace}(f,u0,tspan;
                      callback=nothing,mass_matrix=I) where {isinplace}
    if mass_matrix == I && typeof(f) <: Tuple
      _mm = ((I for i in 1:length(f))...)
    else
      _mm = mass_matrix
    end
    new{typeof(u0),promote_type(map(typeof,tspan)...),
       isinplace,typeof(f),typeof(callback),typeof(_mm)}(
       f,u0,tspan,callback,_mm)
  end
end

function ODEProblem(f,u0,tspan;kwargs...)
  iip = typeof(f)<: Tuple ? isinplace(f[2],3) : isinplace(f,3)
  ODEProblem{iip}(f,u0,tspan;kwargs...)
end

# u'' = f(t,u,du,ddu)
function SecondOrderODEProblem(f,u0,du0,tspan;iip = isinplace(f,4),kwargs...)
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
  ODEProblem{iip}(_f,_u0,tspan;kwargs...)
end
