struct StandardODEProblem end

# Mu' = f
struct ODEProblem{uType,tType,isinplace,P,F,J,C,MM,PT} <:
               AbstractODEProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  p::P
  jac_prototype::J
  callback::C
  mass_matrix::MM
  problem_type::PT
  function ODEProblem{iip}(f,u0,tspan,p=nothing,
                      problem_type=StandardODEProblem();
                      jac_prototype = nothing,
                      callback=nothing,mass_matrix=I) where {iip}
    if mass_matrix == I && typeof(f) <: Tuple
      _mm = ((I for i in 1:length(f))...)
    else
      _mm = mass_matrix
    end
    new{typeof(u0),promote_type(map(typeof,tspan)...),
       iip,typeof(p),typeof(f),typeof(jac_prototype),
       typeof(callback),typeof(_mm),
       typeof(problem_type)}(
       f,u0,tspan,p,jac_prototype,callback,_mm,problem_type)
  end
end

function ODEProblem(f,u0,tspan,p=nothing;kwargs...)
  iip = typeof(f)<: Tuple ? isinplace(f[2],4) : isinplace(f,4)
  ODEProblem{iip}(f,u0,tspan,p;kwargs...)
end

abstract type AbstractDynamicalODEProblem end

struct DynamicalODEProblem{iip} <: AbstractDynamicalODEProblem end
# u' = f1(v)
# v' = f2(t,u)

struct DynamicalODEFunction{iip,F1,F2} <: Function
    f1::F1
    f2::F2
    DynamicalODEFunction{iip}(f1,f2) where iip =
                        new{iip,typeof(f1),typeof(f2)}(f1,f2)
end
function (f::DynamicalODEFunction)(u,p,t)
    ArrayPartition(f.f1(u.x[1],u.x[2],p,t),f.f2(u.x[1],u.x[2],p,t))
end
function (f::DynamicalODEFunction)(du,u,p,t)
    f.f1(du.x[1],u.x[1],u.x[2],p,t)
    f.f2(du.x[2],u.x[1],u.x[2],p,t)
end

function DynamicalODEProblem(f1,f2,du0,u0,tspan,p=nothing;kwargs...)
  iip = isinplace(f1,5)
  DynamicalODEProblem{iip}(f1,f2,du0,u0,tspan,p;kwargs...)
end
function DynamicalODEProblem{iip}(f1,f2,du0,u0,tspan,p=nothing;kwargs...) where iip
    ODEProblem{iip}(DynamicalODEFunction{iip}(f1,f2),(du0,u0),tspan,p;kwargs...)
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
  ODEProblem{iip}(DynamicalODEFunction{iip}(f,f2),_u0,tspan,p,
                  SecondOrderODEProblem{iip}();kwargs...)
end

struct SplitFunction{iip,F1,F2,C} <: Function
    f1::F1
    f2::F2
    cache::C
    SplitFunction{iip}(f1,f2,cache) where iip =
                        new{iip,typeof(f1),typeof(f2),typeof(cache)}(f1,f2,cache)
end
function (f::SplitFunction)(u,p,t)
    f.f1(u,p,t) + f.f2(u,p,t)
end
function (f::SplitFunction)(du,u,p,t)
    f.f1(f.cache,u,p,t)
    f.f2(du,u,p,t)
    du .+= f.cache
end

abstract type AbstractSplitODEProblem end
struct SplitODEProblem{iip} <: AbstractSplitODEProblem end
# u' = Au + f
function SplitODEProblem(f1,f2,u0,tspan,p=nothing;kwargs...)
  iip = isinplace(f2,4)
  SplitODEProblem{iip}(f1,f2,u0,tspan,p;kwargs...)
end
function SplitODEProblem{iip}(f1,f2,u0,tspan,p=nothing;
                                     func_cache=nothing,kwargs...) where iip
  iip ? _func_cache = similar(u0) : _func_cache = nothing
  ODEProblem{iip}(SplitFunction{iip}(f1,f2,_func_cache),u0,tspan,p,
                                     SplitODEProblem{iip}();kwargs...)
end

"""
    set_u0(prob::ODEProblem, u0) -> newprob
Create a new problem from `prob` that has instead initial state `u0`.
"""
function set_u0(prob::ODEProblem, u0)
    ODEProblem{isinplace(prob)}(
    prob.f,
    u0,
    prob.tspan,
    prob.p,
    prob.problem_type;
    jac_prototype = prob.jac_prototype,
    callback = prob.callback,
    mass_matrix = prob.mass_matrix)
end

export set_u0
