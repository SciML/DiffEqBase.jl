const RECOMPILE_BY_DEFAULT = true

abstract type AbstractODEFunction{iip} <: AbstractDiffEqFunction{iip} end
struct ODEFunction{iip,F,Ta,Tt,TJ,JP,TW,TWt,TPJ,S} <: AbstractODEFunction{iip}
  f::F
  analytic::Ta
  tgrad::Tt
  jac::TJ
  jac_prototype::JP
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

struct DynamicalODEFunction{iip,F1,F2,Ta} <: AbstractODEFunction{iip}
  f1::F1
  f2::F2
  analytic::Ta
end

abstract type AbstractDDEFunction{iip} <: AbstractDiffEqFunction{iip} end
struct DDEFunction{iip,F,Ta,Tt,TJ,JP,TW,TWt,TPJ,S} <: AbstractDDEFunction{iip}
  f::F
  analytic::Ta
  tgrad::Tt
  jac::TJ
  jac_prototype::JP
  invW::TW
  invW_t::TWt
  paramjac::TPJ
  syms::S
end

abstract type AbstractDiscreteFunction{iip} <: AbstractDiffEqFunction{iip} end
struct DiscreteFunction{iip,F,Ta} <: AbstractDiscreteFunction{iip}
  f::F
  analytic::Ta
end

abstract type AbstractSDEFunction{iip} <: AbstractDiffEqFunction{iip} end
struct SDEFunction{iip,F,G,Ta,Tt,TJ,JP,TW,TWt,TPJ,S,GG} <: AbstractSDEFunction{iip}
  f::F
  g::G
  analytic::Ta
  tgrad::Tt
  jac::TJ
  jac_prototype::JP
  invW::TW
  invW_t::TWt
  paramjac::TPJ
  ggprime::GG
  syms::S
end

abstract type AbstractRODEFunction{iip} <: AbstractDiffEqFunction{iip} end
struct RODEFunction{iip,F,Ta,Tt,TJ,JP,TW,TWt,TPJ,S} <: AbstractRODEFunction{iip}
  f::F
  analytic::Ta
  tgrad::Tt
  jac::TJ
  jac_prototype::JP
  invW::TW
  invW_t::TWt
  paramjac::TPJ
  syms::S
end

abstract type AbstractDAEFunction{iip} <: AbstractDiffEqFunction{iip} end
struct DAEFunction{iip,F,Ta,Tt,TJ,JP,TW,TWt,TPJ,S} <: AbstractDAEFunction{iip}
  f::F
  analytic::Ta
  tgrad::Tt
  jac::TJ
  jac_prototype::JP
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

function (f::DynamicalODEFunction)(u,p,t)
  ArrayPartition(f.f1(u.x[1],u.x[2],p,t),f.f2(u.x[1],u.x[2],p,t))
end
function (f::DynamicalODEFunction)(du,u,p,t)
  f.f1(du.x[1],u.x[1],u.x[2],p,t)
  f.f2(du.x[2],u.x[1],u.x[2],p,t)
end

(f::SplitFunction)(u,p,t) = f.f1(u,p,t) + f.f2(u,p,t)
(f::SplitFunction)(::Type{Val{:analytic}},args...) = f.analytic(args...)
function (f::SplitFunction)(du,u,p,t)
  f.f1(f.cache,u,p,t)
  f.f2(du,u,p,t)
  du .+= f.cache
end

(f::DiscreteFunction)(args...) = f.f(args...)
(f::DiscreteFunction)(::Type{Val{:analytic}},args...) = f.analytic(args...)

(f::DAEFunction)(args...) = f.f(args...)
(f::DAEFunction)(::Type{Val{:analytic}},args...) = f.analytic(args...)
(f::DAEFunction)(::Type{Val{:tgrad}},args...) = f.tgrad(args...)
(f::DAEFunction)(::Type{Val{:jac}},args...) = f.jac(args...)
(f::DAEFunction)(::Type{Val{:invW}},args...) = f.invW(args...)
(f::DAEFunction)(::Type{Val{:invW_t}},args...) = f.invW_t(args...)
(f::DAEFunction)(::Type{Val{:paramjac}},args...) = f.paramjac(args...)

(f::DDEFunction)(args...) = f.f(args...)
(f::DDEFunction)(::Type{Val{:analytic}},args...) = f.analytic(args...)

(f::SDEFunction)(args...) = f.f(args...)
(f::SDEFunction)(::Type{Val{:analytic}},args...) = f.analytic(args...)
(f::SDEFunction)(::Type{Val{:tgrad}},args...) = f.tgrad(args...)
(f::SDEFunction)(::Type{Val{:jac}},args...) = f.jac(args...)
(f::SDEFunction)(::Type{Val{:invW}},args...) = f.invW(args...)
(f::SDEFunction)(::Type{Val{:invW_t}},args...) = f.invW_t(args...)
(f::SDEFunction)(::Type{Val{:paramjac}},args...) = f.paramjac(args...)

(f::RODEFunction)(args...) = f.f(args...)

######### Basic Constructor

function ODEFunction{iip,true}(f;
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 jac_prototype=nothing,
                 invW=nothing,
                 invW_t=nothing,
                 paramjac = nothing,
                 syms = nothing) where iip
                 if jac == nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                   jac = update_coefficients!
                 end
                 ODEFunction{iip,typeof(f),typeof(analytic),typeof(tgrad),
                 typeof(jac),typeof(jac_prototype),typeof(invW),typeof(invW_t),
                 typeof(paramjac),typeof(syms)}(
                 f,analytic,tgrad,jac,jac_prototype,invW,invW_t,
                 paramjac,syms)
end
function ODEFunction{iip,false}(f;
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 jac_prototype=nothing,
                 invW=nothing,
                 invW_t=nothing,
                 paramjac = nothing,
                 syms = nothing) where iip
                 if jac == nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                   jac = update_coefficients!
                 end
                 ODEFunction{iip,Any,Any,Any,
                 Any,Any,Any,Any,
                 Any,typeof(syms)}(
                 f,analytic,tgrad,jac,jac_prototype,invW,invW_t,
                 paramjac,syms)
end
ODEFunction(f; kwargs...) = ODEFunction{isinplace(f, 4),RECOMPILE_BY_DEFAULT}(f; kwargs...)

@add_kwonly function SplitFunction(f1,f2,cache,analytic)
  f1 = typeof(f1) <: AbstractDiffEqOperator ? f1 : ODEFunction(f1)
  f2 = ODEFunction(f2)
  SplitFunction{isinplace(f2),typeof(f1),typeof(f2),
              typeof(cache),typeof(analytic)}(f1,f2,cache,analytic)
end
SplitFunction{iip,true}(f1,f2; _func_cache=nothing,analytic=nothing) where iip =
SplitFunction{iip,typeof(f1),typeof(f2),typeof(_func_cache),typeof(analytic)}(f1,f2,_func_cache,analytic)
SplitFunction{iip,false}(f1,f2; _func_cache=nothing,analytic=nothing) where iip =
SplitFunction{iip,Any,Any,Any}(f1,f2,_func_cache,analytic)
SplitFunction(f1,f2; kwargs...) = SplitFunction{isinplace(f2, 4)}(f1, f2; kwargs...)
SplitFunction{iip}(f1,f2; kwargs...) where iip =
SplitFunction{iip,RECOMPILE_BY_DEFAULT}(
typeof(f1) <: AbstractDiffEqOperator ? f1 : ODEFunction(f1),
ODEFunction{iip}(f2); kwargs...)

@add_kwonly function DynamicalODEFunction{iip}(f1,f2,analytic) where iip
  f1 = ODEFunction(f1)
  f2 != nothing && (f2 = ODEFunction(f2))
  DynamicalODEFunction{iip,typeof(f1),typeof(f2),typeof(analytic)}(f1,f2,analytic)
end
DynamicalODEFunction{iip,true}(f1,f2;analytic=nothing) where iip =
DynamicalODEFunction{iip,typeof(f1),typeof(f2),typeof(analytic)}(f1,f2,analytic)
DynamicalODEFunction{iip,false}(f1,f2;analytic=nothing) where iip =
DynamicalODEFunction{iip,Any,Any,Any}(f1,f2,analytic)
DynamicalODEFunction(f1,f2=nothing; kwargs...) = DynamicalODEFunction{isinplace(f1, 5)}(f1, f2; kwargs...)
DynamicalODEFunction{iip}(f1,f2; kwargs...) where iip =
DynamicalODEFunction{iip,RECOMPILE_BY_DEFAULT}(ODEFunction{iip}(f1), ODEFunction{iip}(f2); kwargs...)

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

function SDEFunction{iip,true}(f,g;
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 jac_prototype=nothing,
                 invW=nothing,
                 invW_t=nothing,
                 paramjac = nothing,
                 ggprime = nothing,
                 syms = nothing) where iip
                 if jac == nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                   jac = update_coefficients!
                 end
                 SDEFunction{iip,typeof(f),typeof(g),
                 typeof(analytic),typeof(tgrad),
                 typeof(jac),typeof(jac_prototype),typeof(invW),typeof(invW_t),
                 typeof(paramjac),typeof(syms),
                 typeof(ggprime)}(
                 f,g,analytic,tgrad,jac,jac_prototype,invW,invW_t,
                 paramjac,ggprime,syms)
end
function SDEFunction{iip,false}(f,g;
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 jac_prototype=nothing,
                 invW=nothing,
                 invW_t=nothing,
                 paramjac = nothing,
                 ggprime = nothing,
                 syms = nothing) where iip
                 if jac == nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                   jac = update_coefficients!
                 end
                 SDEFunction{iip,Any,Any,Any,Any,
                 Any,Any,Any,Any,
                 Any,typeof(syms),Any}(
                 f,g,analytic,tgrad,jac,jac_prototype,invW,invW_t,
                 paramjac,ggprime,syms)
end
SDEFunction(f,g; kwargs...) = SDEFunction{isinplace(f, 4),RECOMPILE_BY_DEFAULT}(f,g; kwargs...)

function RODEFunction{iip,true}(f;
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 jac_prototype=nothing,
                 invW=nothing,
                 invW_t=nothing,
                 paramjac = nothing,
                 syms = nothing) where iip
                 if jac == nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                   jac = update_coefficients!
                 end
                 RODEFunction{iip,typeof(f),
                 typeof(analytic),typeof(tgrad),
                 typeof(jac),typeof(jac_prototype),typeof(invW),typeof(invW_t),
                 typeof(paramjac),typeof(syms)}(
                 f,analytic,tgrad,jac,jac_prototype,invW,invW_t,
                 paramjac,syms)
end
function RODEFunction{iip,false}(f;
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 invW=nothing,
                 invW_t=nothing,
                 paramjac = nothing,
                 syms = nothing) where iip
                 if jac == nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                   jac = update_coefficients!
                 end
                 RODEFunction{iip,Any,Any,Any,
                 Any,Any,Any,Any,
                 Any,typeof(syms)}(
                 f,analytic,tgrad,jac,jac_prototype,invW,invW_t,
                 paramjac,syms)
end
RODEFunction(f; kwargs...) = RODEFunction{isinplace(f, 5),RECOMPILE_BY_DEFAULT}(f; kwargs...)

function DAEFunction{iip,true}(f;
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 jac_prototype=nothing,
                 invW=nothing,
                 invW_t=nothing,
                 paramjac = nothing,
                 syms = nothing) where iip
                 if jac == nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                   jac = update_coefficients!
                 end
                 DAEFunction{iip,typeof(f),typeof(analytic),typeof(tgrad),
                 typeof(jac),typeof(jac_prototype),typeof(invW),typeof(invW_t),
                 typeof(paramjac),typeof(syms)}(
                 f,analytic,tgrad,jac,jac_prototype,invW,invW_t,
                 paramjac,syms)
end
function DAEFunction{iip,false}(f;
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 jac_prototype=nothing,
                 invW=nothing,
                 invW_t=nothing,
                 paramjac = nothing,
                 syms = nothing) where iip
                 if jac == nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                   jac = update_coefficients!
                 end
                 DAEFunction{iip,Any,Any,Any,
                 Any,Any,Any,Any,
                 Any,typeof(syms)}(
                 f,analytic,tgrad,jac,jac_prototype,invW,invW_t,
                 paramjac,syms)
end
DAEFunction(f; kwargs...) = DAEFunction{isinplace(f, 5),RECOMPILE_BY_DEFAULT}(f; kwargs...)

function DDEFunction{iip,true}(f;
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 jac_prototype=nothing,
                 invW=nothing,
                 invW_t=nothing,
                 paramjac = nothing,
                 syms = nothing) where iip
                 if jac == nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                   jac = update_coefficients!
                 end
                 DDEFunction{iip,typeof(f),typeof(analytic),typeof(tgrad),
                 typeof(jac),typeof(jac_prototype),typeof(invW),typeof(invW_t),
                 typeof(paramjac),typeof(syms)}(
                 f,analytic,tgrad,jac,jac_prototype,invW,invW_t,
                 paramjac,syms)
end
function DDEFunction{iip,false}(f;
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 jac_prototype=nothing,
                 invW=nothing,
                 invW_t=nothing,
                 paramjac = nothing,
                 syms = nothing) where iip
                 if jac == nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                   jac = update_coefficients!
                 end
                 DDEFunction{iip,Any,Any,Any,
                 Any,Any,Any,Any,
                 Any,typeof(syms)}(
                 f,analytic,tgrad,jac,jac_prototype,invW,invW_t,
                 paramjac,syms)
end
DDEFunction(f; kwargs...) = DDEFunction{isinplace(f, 5),RECOMPILE_BY_DEFAULT}(f; kwargs...)

########## Existance Functions

has_analytic(f::AbstractDiffEqFunction) = f.analytic != nothing
has_jac(f::AbstractDiffEqFunction) = f.jac != nothing
has_tgrad(f::AbstractDiffEqFunction) = f.tgrad != nothing
has_invW(f::AbstractDiffEqFunction) = f.invW != nothing
has_invW_t(f::AbstractDiffEqFunction) = f.invW_t != nothing
has_paramjac(f::AbstractDiffEqFunction) = f.paramjac != nothing
has_syms(f::AbstractDiffEqFunction) = f.syms != nothing

# TODO: find an appropriate way to check `has_*`
has_jac(f::SplitFunction) = f.f1.jac != nothing
has_tgrad(f::SplitFunction) = f.f1.tgrad != nothing
has_invW(f::SplitFunction) = f.f1.invW != nothing
has_invW_t(f::SplitFunction) = f.f1.invW_t != nothing
has_paramjac(f::SplitFunction) = f.f1.paramjac != nothing

has_jac(f::DynamicalODEFunction) = f.f1.jac != nothing
has_tgrad(f::DynamicalODEFunction) = f.f1.tgrad != nothing
has_invW(f::DynamicalODEFunction) = f.f1.invW != nothing
has_invW_t(f::DynamicalODEFunction) = f.f1.invW_t != nothing
has_paramjac(f::DynamicalODEFunction) = f.f1.paramjac != nothing

######### Compatibility Constructor from Tratis

ODEFunction{iip}(f::T) where {iip,T} = return T<:ODEFunction ? f : convert(ODEFunction{iip},f)
function Base.convert(::Type{ODEFunction},f)
  if __has_analytic(f)
    analytic = (args...) -> f(Val{:analytic},args...)
  else
    analytic = nothing
  end
  if __has_jac(f)
    @warn("The overloading form for Jacobians is deprecated. Use the DiffEqFunction")
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
  ODEFunction(f;analytic=analytic,tgrad=tgrad,jac=jac,invW=invW,
              invW_t=invW_t,paramjac=paramjac,syms=syms)
end
function Base.convert(::Type{ODEFunction{iip}},f) where iip
  if __has_analytic(f)
    analytic = (args...) -> f(Val{:analytic},args...)
  else
    analytic = nothing
  end
  if __has_jac(f)
    @warn("The overloading form for Jacobians is deprecated. Use the DiffEqFunction")
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
  ODEFunction{iip,RECOMPILE_BY_DEFAULT}(f;analytic=analytic,tgrad=tgrad,jac=jac,invW=invW,
              invW_t=invW_t,paramjac=paramjac,syms=syms)
end

DiscreteFunction{iip}(f::T) where {iip,T} = return T<:DiscreteFunction ? f : convert(DiscreteFunction{iip},f)
function Base.convert(::Type{DiscreteFunction},f)
  if __has_analytic(f)
    analytic = (args...) -> f(Val{:analytic},args...)
  else
    analytic = nothing
  end
  DiscreteFunction(f;analytic=analytic)
end
function Base.convert(::Type{DiscreteFunction{iip}},f) where iip
  if __has_analytic(f)
    analytic = (args...) -> f(Val{:analytic},args...)
  else
    analytic = nothing
  end
  DiscreteFunction{iip,RECOMPILE_BY_DEFAULT}(f;analytic=analytic)
end

DAEFunction{iip}(f::T) where {iip,T} = return T<:DAEFunction ? f : convert(DAEFunction{iip},f)
function Base.convert(::Type{DAEFunction},f)
  if __has_analytic(f)
    analytic = (args...) -> f(Val{:analytic},args...)
  else
    analytic = nothing
  end
  if __has_jac(f)
    @warn("The overloading form for Jacobians is deprecated. Use the DiffEqFunction")
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
  DAEFunction(f;analytic=analytic,tgrad=tgrad,jac=jac,invW=invW,
              invW_t=invW_t,paramjac=paramjac,syms=syms)
end
function Base.convert(::Type{DAEFunction{iip}},f) where iip
  if __has_analytic(f)
    analytic = (args...) -> f(Val{:analytic},args...)
  else
    analytic = nothing
  end
  if __has_jac(f)
    @warn("The overloading form for Jacobians is deprecated. Use the DiffEqFunction")
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
  DAEFunction{iip,RECOMPILE_BY_DEFAULT}(f;analytic=analytic,tgrad=tgrad,jac=jac,invW=invW,
              invW_t=invW_t,paramjac=paramjac,syms=syms)
end

DDEFunction{iip}(f::T) where {iip,T} = return T<:DDEFunction ? f : convert(DDEFunction{iip},f)
function Base.convert(::Type{DDEFunction},f)
  if __has_analytic(f)
    analytic = (args...) -> f(Val{:analytic},args...)
  else
    analytic = nothing
  end
  if __has_syms(f)
    syms = f.syms
  else
    syms = nothing
  end
  DDEFunction(f;analytic=analytic,syms=syms)
end
function Base.convert(::Type{DDEFunction{iip}},f) where iip
  if __has_analytic(f)
    analytic = (args...) -> f(Val{:analytic},args...)
  else
    analytic = nothing
  end
  if __has_syms(f)
    syms = f.syms
  else
    syms = nothing
  end
  DDEFunction{iip,RECOMPILE_BY_DEFAULT}(f;analytic=analytic,syms=syms)
end

SDEFunction{iip}(f::T,g::T2) where {iip,T,T2} = return T<:SDEFunction ? f : convert(SDEFunction{iip},f,g)
function Base.convert(::Type{SDEFunction},f,g)
  if __has_analytic(f)
    analytic = (args...) -> f(Val{:analytic},args...)
  else
    analytic = nothing
  end
  if __has_jac(f)
    @warn("The overloading form for Jacobians is deprecated. Use the DiffEqFunction")
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
  SDEFunction(f,g;analytic=analytic,tgrad=tgrad,jac=jac,invW=invW,
              invW_t=invW_t,paramjac=paramjac,syms=syms)
end
function Base.convert(::Type{SDEFunction{iip}},f,g) where iip
  if __has_analytic(f)
    analytic = (args...) -> f(Val{:analytic},args...)
  else
    analytic = nothing
  end
  if __has_jac(f)
    @warn("The overloading form for Jacobians is deprecated. Use the DiffEqFunction")
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
  SDEFunction{iip,RECOMPILE_BY_DEFAULT}(f,g;analytic=analytic,
              tgrad=tgrad,jac=jac,invW=invW,
              invW_t=invW_t,paramjac=paramjac,syms=syms)
end

RODEFunction{iip}(f::T) where {iip,T,T2} = return T<:RODEFunction ? f : convert(RODEFunction{iip},f)
function Base.convert(::Type{RODEFunction},f)
  if __has_analytic(f)
    analytic = (args...) -> f(Val{:analytic},args...)
  else
    analytic = nothing
  end
  if __has_jac(f)
    @warn("The overloading form for Jacobians is deprecated. Use the DiffEqFunction")
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
  RODEFunction(f;analytic=analytic,tgrad=tgrad,jac=jac,invW=invW,
              invW_t=invW_t,paramjac=paramjac,syms=syms)
end
function Base.convert(::Type{RODEFunction{iip}},f) where iip
  if __has_analytic(f)
    analytic = (args...) -> f(Val{:analytic},args...)
  else
    analytic = nothing
  end
  if __has_jac(f)
    @warn("The overloading form for Jacobians is deprecated. Use the DiffEqFunction")
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
  RODEFunction{iip,RECOMPILE_BY_DEFAULT}(f;analytic=analytic,
              tgrad=tgrad,jac=jac,invW=invW,
              invW_t=invW_t,paramjac=paramjac,syms=syms)
end
