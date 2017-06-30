# Mu' = f
type ODEProblem{uType,tType,isinplace,F,C,MM} <: AbstractODEProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
  mass_matrix::MM
end

function ODEProblem(f,u0,tspan;
                    iip = typeof(f)<: Tuple ? isinplace(f[2],3) : isinplace(f,3),callback=nothing,mass_matrix=I)
  if mass_matrix == I && typeof(f) <: Tuple
    _mm = ((I for i in 1:length(f))...)
  else
    _mm = mass_matrix
  end
  ODEProblem{typeof(u0),promote_type(map(typeof,tspan)...),
             iip,typeof(f),typeof(callback),typeof(_mm)}(
             f,u0,tspan,callback,_mm)
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
  ODEProblem{typeof(_u0),promote_type(map(typeof,tspan)...),
             iip,typeof(_f),typeof(callback),typeof(_mass_matrix)}(
             _f,_u0,tspan,callback,_mass_matrix)
end
