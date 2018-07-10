const RECOMPILE_BY_DEFAULT = true

abstract type AbstractODEFunction{iip} <: AbstractDiffEqFunction{iip} end
struct ODEFunction{iip,F,Ta,Tt,TJ,TW,TWt,TPJ,S} <: AbstractODEFunction{iip}
  f::F
  analytic::Ta
  tgrad::Tt
  jac::TJ
  invW::TW
  invW_t::TWt
  paramjac::TPJ
  syms::S
end

struct SplitFunction{iip,F1,F2,C,Ta} <: AbstractODEFunction{iip}
    f1::F1
    f2::F2
    cache::C
    analytic::Ta
end

abstract type AbstractDiscreteFunction{iip} <: AbstractDiffEqFunction{iip} end
struct DiscreteFunction{iip,F,Ta} <: AbstractDiscreteFunction{iip}
  f::F
  analytic::Ta
end

######### Backwards Compatibility Overloads

(f::ODEFunction)(args...) = f.f(args...)
(f::ODEFunction)(::Type{Val{:analytic}},args...) = f.analytic(args...)
(f::ODEFunction)(::Type{Val{:tgrad}},args...) = f.tgrad(args...)
(f::ODEFunction)(::Type{Val{:jac}},args...) = f.jac(args...)
(f::ODEFunction)(::Type{Val{:invW}},args...) = f.invW(args...)
(f::ODEFunction)(::Type{Val{:invW_t}},args...) = f.invW_t(args...)
(f::ODEFunction)(::Type{Val{:paramjac}},args...) = f.paramjac(args...)

(f::SplitFunction)(u,p,t) = f.f1(u,p,t) + f.f2(u,p,t)
(f::SplitFunction)(::Type{Val{:analytic}},args...) = f.analytic(args...)
function (f::SplitFunction)(du,u,p,t)
    f.f1(f.cache,u,p,t)
    f.f2(du,u,p,t)
    du .+= f.cache
end

(f::DiscreteFunction)(args...) = f.f(args...)
(f::DiscreteFunction)(::Type{Val{:analytic}},args...) = f.analytic(args...)

######### Basic Constructor

function ODEFunction{iip,true}(f;
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 invW=nothing,
                 invW_t=nothing,
                 paramjac = nothing,
                 syms = nothing) where iip
                 ODEFunction{iip,typeof(f),typeof(analytic),typeof(tgrad),
                 typeof(jac),typeof(invW),typeof(invW_t),
                 typeof(paramjac),typeof(syms)}(
                 f,analytic,tgrad,jac,invW,invW_t,
                 paramjac,syms)
end
function ODEFunction{iip,false}(f;
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 invW=nothing,
                 invW_t=nothing,
                 paramjac = nothing,
                 syms = nothing) where iip
                 ODEFunction{iip,Any,Any,Any,
                 Any,Any,Any,
                 Any,typeof(syms)}(
                 f,analytic,tgrad,jac,invW,invW_t,
                 paramjac,syms)
end
ODEFunction(f; kwargs...) = ODEFunction{isinplace(f, 4),RECOMPILE_BY_DEFAULT}(f; kwargs...)

SplitFunction{iip,true}(f1,f2; _func_cache=nothing,analytic=nothing) where iip =
  SplitFunction{iip,typeof(f1),typeof(f2),typeof(_func_cache),typeof(analytic)}(f1,f2,_func_cache,analytic)
SplitFunction{iip,false}(f1,f2; _func_cache=nothing,analytic=nothing) where iip =
  SplitFunction{iip,Any,Any,Any}(f1,f2,_func_cache,analytic)
function SplitFunction(f1,f2; kwargs...)
  iip2 = isinplace(f2, 4)
  SplitFunction{iip2,RECOMPILE_BY_DEFAULT}(f1, f2; kwargs...)
end

function DiscreteFunction{iip,true}(f;
                 analytic=nothing) where iip
                 DiscreteFunction{iip,typeof(f),typeof(analytic)}(
                 f,analytic)
end
function DiscreteFunction{iip,false}(f;
                 analytic=nothing) where iip
                 DiscreteFunction{iip,Any,Any}(
                 f,analytic)
end
DiscreteFunction(f; kwargs...) = DiscreteFunction{isinplace(f, 4),RECOMPILE_BY_DEFAULT}(f; kwargs...)

########## Existance Functions

has_jac(f::ODEFunction) = f.jac != nothing
has_analytic(f::ODEFunction) = f.analytic != nothing
has_tgrad(f::ODEFunction) = f.tgrad != nothing
has_invW(f::ODEFunction) = f.invW != nothing
has_invW_t(f::ODEFunction) = f.invW_t != nothing
has_paramjac(f::ODEFunction) = f.paramjac != nothing
has_syms(f::ODEFunction) = f.syms != nothing

has_analytic(f::SplitFunction) = f.analytic != nothing

has_analytic(f::DiscreteFunction) = f.analytic != nothing

######### Compatibility Constructor from Tratis

function Base.convert(::Type{ODEFunction},f)
  if __has_analytic(f)
    analytic = (args...) -> f(Val{:analytic},args...)
  else
    analytic = nothing
  end
  if __has_jac(f)
    warn("The overloading form for Jacobians is deprecated. Use the DiffEqFunction")
    jac = (args...) -> f(Val{:jac},args...)
  else
    jac = nothing
  end
  if __has_tgrad(f)
    tgrad = (args...) -> f(Val{:tgrad},args...)
  else
    tgrad = nothing
  end
  if __has_invW(f)
    invW = (args...) -> f(Val{:invW},args...)
  else
    invW = nothing
  end
  if __has_invW_t(f)
    invW_t = (args...) -> f(Val{:invW_t},args...)
  else
    invW_t = nothing
  end
  if __has_paramjac(f)
    paramjac = (args...) -> f(Val{:paramjac},args...)
  else
    paramjac = nothing
  end
  if __has_syms(f)
    syms = f.syms
  else
    syms = nothing
  end
  ODEFunction(f,analytic=analytic,tgrad=tgrad,jac=jac,invW=invW,
              invW_t=invW_t,paramjac=paramjac,syms=syms)
end
function Base.convert(::Type{ODEFunction{iip}},f) where iip
  if __has_analytic(f)
    analytic = (args...) -> f(Val{:analytic},args...)
  else
    analytic = nothing
  end
  if __has_jac(f)
    warn("The overloading form for Jacobians is deprecated. Use the DiffEqFunction")
    jac = (args...) -> f(Val{:jac},args...)
  else
    jac = nothing
  end
  if __has_tgrad(f)
    tgrad = (args...) -> f(Val{:tgrad},args...)
  else
    tgrad = nothing
  end
  if __has_invW(f)
    invW = (args...) -> f(Val{:invW},args...)
  else
    invW = nothing
  end
  if __has_invW_t(f)
    invW_t = (args...) -> f(Val{:invW_t},args...)
  else
    invW_t = nothing
  end
  if __has_paramjac(f)
    paramjac = (args...) -> f(Val{:paramjac},args...)
  else
    paramjac = nothing
  end
  if __has_syms(f)
    syms = f.syms
  else
    syms = nothing
  end
  ODEFunction{iip,RECOMPILE_BY_DEFAULT}(f,analytic=analytic,tgrad=tgrad,jac=jac,invW=invW,
              invW_t=invW_t,paramjac=paramjac,syms=syms)
end

function Base.convert(::Type{DiscreteFunction},f)
  if __has_analytic(f)
    analytic = (args...) -> f(Val{:analytic},args...)
  else
    analytic = nothing
  end
  DiscreteFunction(f,analytic=analytic)
end
function Base.convert(::Type{DiscreteFunction{iip}},f) where iip
  if __has_analytic(f)
    analytic = (args...) -> f(Val{:analytic},args...)
  else
    analytic = nothing
  end
  DiscreteFunction{iip,RECOMPILE_BY_DEFAULT}(f,analytic=analytic)
end
