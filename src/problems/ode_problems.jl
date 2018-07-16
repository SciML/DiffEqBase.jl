struct StandardODEProblem end

# Mu' = f
struct ODEProblem{uType,tType,isinplace,P,F,J,C,MM,PT} <:
               AbstractODEProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::tType
  p::P
  jac_prototype::J
  callback::C
  mass_matrix::MM
  problem_type::PT
  @add_kwonly function ODEProblem(f::AbstractODEFunction,u0,tspan,p=nothing,
                      problem_type=StandardODEProblem();
                      jac_prototype = nothing,
                      callback=nothing,mass_matrix=I)
    _tspan = promote_tspan(tspan)
    if mass_matrix == I && typeof(f) <: Tuple
      _mm = ((I for i in 1:length(f))...,)
    else
      _mm = mass_matrix
    end
    new{typeof(u0),typeof(_tspan),
       isinplace(f),typeof(p),typeof(f),typeof(jac_prototype),
       typeof(callback),typeof(_mm),
       typeof(problem_type)}(
       f,u0,_tspan,p,jac_prototype,callback,_mm,problem_type)
  end

  function ODEProblem{iip}(f,u0,tspan,p=nothing;kwargs...) where {iip}
    ODEProblem(convert(ODEFunction{iip},f),u0,tspan,p;kwargs...)
  end
end

function ODEProblem(f,u0,tspan,p=nothing;kwargs...)
  ODEProblem(convert(ODEFunction,f),u0,tspan,p;kwargs...)
end

abstract type AbstractDynamicalODEProblem end

struct DynamicalODEProblem{iip} <: AbstractDynamicalODEProblem end
# u' = f1(v)
# v' = f2(t,u)
function DynamicalODEProblem(f::DynamicalODEFunction,du0,u0,tspan,p=nothing;mass_matrix=(I,I),kwargs...)
  ODEProblem(f,(du0,u0),tspan,p;mass_matrix=mass_matrix,kwargs...)
end
function DynamicalODEProblem(f1,f2,du0,u0,tspan,p=nothing;mass_matrix=(I,I),kwargs...)
  ODEProblem(DynamicalODEFunction(f1,f2),(du0,u0),tspan,p;mass_matrix=mass_matrix,kwargs...)
end
function DynamicalODEProblem{iip}(f1,f2,du0,u0,tspan,p=nothing;mass_matrix=(I,I),kwargs...) where iip
  ODEProblem(DynamicalODEFunction{iip}(f1,f2),(du0,u0),tspan,p;mass_matrix=mass_matrix,kwargs...)
end

# u'' = f(t,u,du,ddu)
struct SecondOrderODEProblem{iip} <: AbstractDynamicalODEProblem end
function SecondOrderODEProblem(f,du0,u0,tspan,p=nothing;kwargs...)
  iip = isinplace(f,5)
  SecondOrderODEProblem{iip}(f,du0,u0,tspan,p;kwargs...)
end
function SecondOrderODEProblem{iip}(f,du0,u0,tspan,p=nothing;kwargs...) where iip
  if iip
    f2 = function (du,v,u,p,t)
      du .= v
    end
  else
    f2 = function (v,u,p,t)
      v
    end
  end
  _u0 = (du0,u0)
  ODEProblem(DynamicalODEFunction{iip}(f,f2),_u0,tspan,p,
                  SecondOrderODEProblem{iip}();kwargs...)
end
function SecondOrderODEProblem(f::DynamicalODEFunction,du0,u0,tspan,p=nothing;kwargs...)
  iip = isinplace(f.f1, 5)
  _u0 = (du0,u0)
  if f.f2.f == nothing
    if iip
      f2 = function (du,v,u,p,t)
        du .= v
      end
    else
      f2 = function (v,u,p,t)
        v
      end
    end
    return ODEProblem(DynamicalODEFunction{iip}(f.f1,f2;analytic=f.analytic),_u0,tspan,p,
                  SecondOrderODEProblem{iip}();kwargs...)
  else
    return ODEProblem(DynamicalODEFunction{iip}(f.f1,f.f2;analytic=f.analytic),_u0,tspan,p,
                  SecondOrderODEProblem{iip}();kwargs...)
  end
end

abstract type AbstractSplitODEProblem end
struct SplitODEProblem{iip} <: AbstractSplitODEProblem end
# u' = Au + f
function SplitODEProblem(f1,f2,u0,tspan,p=nothing;kwargs...)
  f = SplitFunction(f1,f2)
  SplitODEProblem(f,u0,tspan,p;kwargs...)
end
function SplitODEProblem{iip}(f1,f2,u0,tspan,p=nothing;kwargs...) where iip
  f = SplitFunction{iip}(f1,f2)
  SplitODEProblem(f,u0,tspan,p;kwargs...)
end
SplitODEProblem(f::SplitFunction,u0,tspan,p=nothing;kwargs...) =
  SplitODEProblem{isinplace(f)}(f,u0,tspan,p;kwargs...)
function SplitODEProblem{iip}(f::SplitFunction,u0,tspan,p=nothing;kwargs...) where iip
  if f.cache == nothing && iip
    cache = similar(u0)
    f = SplitFunction{iip}(f.f1, f.f2;
                     _func_cache=cache, analytic=f.analytic)
  end
  ODEProblem(f,u0,tspan,p,SplitODEProblem{iip}();kwargs...)
end
