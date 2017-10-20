struct DiffEqFunction{iip,F,Ta,Tt,TJ,TIJ,TW,TWt,TPJ} <: AbstractDiffEqFunction{iip}
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

(f::DiffEqFunction)(args...) = f.f(args...)
(f::DiffEqFunction)(::Type{Val{:analytic}},args...) = f.analytic(args...)
(f::DiffEqFunction)(::Type{Val{:tgrad}},args...) = f.tgrad(args...)
(f::DiffEqFunction)(::Type{Val{:jac}},args...) = f.jac(args...)
(f::DiffEqFunction)(::Type{Val{:invjac}},args...) = f.invjac(args...)
(f::DiffEqFunction)(::Type{Val{:invW}},args...) = f.invW(args...)
(f::DiffEqFunction)(::Type{Val{:invW_t}},args...) = f.invW_t(args...)
(f::DiffEqFunction)(::Type{Val{:paramjac}},args...) = f.paramjac(args...)

######### Basic Constructor

function DiffEqFunction{iip}(f;analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 invjac=nothing,
                 invW=nothing,
                 invW_t=nothing,
                 paramjac = nothing) where iip
                 DiffEqFunction{iip,typeof(f),typeof(analytic),typeof(tgrad),
                 typeof(jac),typeof(invjac),typeof(invW),typeof(invW_t),
                 typeof(paramjac)}(f,analytic,tgrad,jac,invjac,invW,invW_t,
                 paramjac)
end

######### No Specialization Constructors

struct NSODEFunction{iip} <: AbstractDiffEqFunction{iip} end

function NSODEFunction(f,t,u;kwargs...)
  iip = typeof(f)<: Tuple ? isinplace(f[2],3) : isinplace(f,3)
  NSODEFunction{iip}(f,t,u;kwargs...)
end

function NSODEFunction{iip}(f,t,u;analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 invjac=nothing,
                 invW=nothing,
                 invW_t=nothing,
                 paramjac = nothing) where iip
                 if iip
                   _f = (t,u,du) -> (f(t,u,du);nothing)
                   wrap_f = FunctionWrappers.FunctionWrapper{Void,Tuple{typeof(t),typeof(u),typeof(u)}}(_f)
                 else
                   _f = f
                   wrap_f = FunctionWrappers.FunctionWrapper{typeof(u),Tuple{typeof(t),typeof(u)}}(_f)
                 end
                 DiffEqFunction{iip,typeof(wrap_f),typeof(analytic),typeof(tgrad),
                 typeof(jac),typeof(invjac),typeof(invW),typeof(invW_t),
                 typeof(paramjac)}(wrap_f,analytic,tgrad,jac,invjac,invW,invW_t,
                 paramjac)
end

########## Existance Functions

has_jac(f::DiffEqFunction) = f.jac != nothing
has_invjac(f::DiffEqFunction) = f.invjac != nothing
has_analytic(f::DiffEqFunction) = f.analytic != nothing
has_tgrad(f::DiffEqFunction) = f.tgrad != nothing
has_invW(f::DiffEqFunction) = f.invW != nothing
has_invW_t(f::DiffEqFunction) = f.invW_t != nothing
has_paramjac(f::DiffEqFunction) = f.paramjac != nothing

######### Compatibility Constructor from Tratis

function CompatibilityDiffEqFunction(f)
  if has_analytic(f)
    analytic = (args...) -> f(Val{:analytic},args...)
  else
    analytic = nothing
  end
  if has_jac(f)
    jac = (args...) -> f(Val{:jac},args...)
  else
    jac = nothing
  end
  if has_invjac(f)
    invjac = (args...) -> f(Val{:invjac},args...)
  else
    invjac = nothing
  end
  if has_tgrad(f)
    tgrad = (args...) -> f(Val{:tgrad},args...)
  else
    tgrad = nothing
  end
  if has_invW(f)
    invW = (args...) -> f(Val{:invW},args...)
  else
    invW = nothing
  end
  if has_invW_t(f)
    invW_t = (args...) -> f(Val{:invW_t},args...)
  else
    invW_t = nothing
  end
  if has_paramjac(f)
    paramjac = (args...) -> f(Val{:paramjac},args...)
  else
    paramjac = nothing
  end
  DiffEqFunction(f,analytic,tgrad,jac,invjac,invW,invW_t,paramjac)
end
