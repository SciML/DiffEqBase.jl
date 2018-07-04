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

######### Backwards Compatibility Overloads

(f::ODEFunction)(args...) = f.f(args...)
(f::ODEFunction)(::Type{Val{:analytic}},args...) = f.analytic(args...)
(f::ODEFunction)(::Type{Val{:tgrad}},args...) = f.tgrad(args...)
(f::ODEFunction)(::Type{Val{:jac}},args...) = f.jac(args...)
(f::ODEFunction)(::Type{Val{:invW}},args...) = f.invW(args...)
(f::ODEFunction)(::Type{Val{:invW_t}},args...) = f.invW_t(args...)
(f::ODEFunction)(::Type{Val{:paramjac}},args...) = f.paramjac(args...)

######### Basic Constructor

function ODEFunction{iip}(f;
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
ODEFunction(f; kwargs...) = ODEFunction{isinplace(f, 4)}(f; kwargs...)

########## Existance Functions

has_jac(f::ODEFunction) = f.jac != nothing
has_analytic(f::ODEFunction) = f.analytic != nothing
has_tgrad(f::ODEFunction) = f.tgrad != nothing
has_invW(f::ODEFunction) = f.invW != nothing
has_invW_t(f::ODEFunction) = f.invW_t != nothing
has_paramjac(f::ODEFunction) = f.paramjac != nothing
has_syms(f::ODEFunction) = f.syms != nothing

######### Compatibility Constructor from Tratis

function Base.convert(::Type{ODEFunction},f)
  if __has_analytic(f)
    analytic = (args...) -> f(Val{:analytic},args...)
  else
    analytic = nothing
  end
  if __has_jac(f)
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
  ODEFunction{iip}(f,analytic=analytic,tgrad=tgrad,jac=jac,invW=invW,
              invW_t=invW_t,paramjac=paramjac,syms=syms)
end
