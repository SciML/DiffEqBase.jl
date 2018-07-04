struct ODEFunction{iip,F,Ta,Tt,TJ,TIJ,TW,TWt,TPJ} <: AbstractDiffEqFunction{iip}
  f::F
  analytic::Ta
  tgrad::Tt
  jac::TJ
  invjac::TIJ
  invW::TW
  invW_t::TWt
  paramjac::TPJ
end

######### Backwards Compatibility Overloads

(f::ODEFunction)(args...) = f.f(args...)
(f::ODEFunction)(::Type{Val{:analytic}},args...) = f.analytic(args...)
(f::ODEFunction)(::Type{Val{:tgrad}},args...) = f.tgrad(args...)
(f::ODEFunction)(::Type{Val{:jac}},args...) = f.jac(args...)
(f::ODEFunction)(::Type{Val{:invjac}},args...) = f.invjac(args...)
(f::ODEFunction)(::Type{Val{:invW}},args...) = f.invW(args...)
(f::ODEFunction)(::Type{Val{:invW_t}},args...) = f.invW_t(args...)
(f::ODEFunction)(::Type{Val{:paramjac}},args...) = f.paramjac(args...)

######### Basic Constructor

function ODEFunction{iip}(f;
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 invjac=nothing,
                 invW=nothing,
                 invW_t=nothing,
                 paramjac = nothing) where iip
                 ODEFunction{iip,typeof(f),typeof(analytic),typeof(tgrad),
                 typeof(jac),typeof(invjac),typeof(invW),typeof(invW_t),
                 typeof(paramjac)}(f,analytic,tgrad,jac,invjac,invW,invW_t,
                 paramjac)
end
ODEFunction(f; kwargs...) = ODEFunction{isinplace(f, 4)}(f; kwargs...)

########## Existance Functions

has_jac(f::ODEFunction) = f.jac != nothing
has_invjac(f::ODEFunction) = f.invjac != nothing
has_analytic(f::ODEFunction) = f.analytic != nothing
has_tgrad(f::ODEFunction) = f.tgrad != nothing
has_invW(f::ODEFunction) = f.invW != nothing
has_invW_t(f::ODEFunction) = f.invW_t != nothing
has_paramjac(f::ODEFunction) = f.paramjac != nothing

######### Compatibility Constructor from Tratis

function CompatibilityODEFunction(f)
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
  if __has_invjac(f)
    invjac = (args...) -> f(Val{:invjac},args...)
  else
    invjac = nothing
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
  ODEFunction(f,analytic,tgrad,jac,invjac,invW,invW_t,paramjac)
end
