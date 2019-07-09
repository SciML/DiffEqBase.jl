const RECOMPILE_BY_DEFAULT = true
const UNPACK_BY_DEFAULT = true

"""
$(TYPEDEF)

TODO
"""
abstract type AbstractODEFunction{iip} <: AbstractDiffEqFunction{iip} end

"""
$(TYPEDEF)

TODO
"""
struct ODEFunction{iip,unpack,F,TMM,Ta,Tt,TJ,JP,TW,TWt,TPJ,S,TCV} <: AbstractODEFunction{iip}
  f::F
  mass_matrix::TMM
  analytic::Ta
  tgrad::Tt
  jac::TJ
  jac_prototype::JP
  Wfact::TW
  Wfact_t::TWt
  paramjac::TPJ
  syms::S
  colorvec::TCV

  function ODEFunction{iip,recompile,unpack}(f;
                                             mass_matrix = I,
                                             analytic = nothing,
                                             tgrad = nothing,
                                             jac = nothing,
                                             jac_prototype = nothing,
                                             Wfact = nothing,
                                             Wfact_t = nothing,
                                             paramjac = nothing,
                                             syms = nothing,
                                             colorvec = nothing) where {iip,recompile,unpack}
    if mass_matrix === I && typeof(f) <: Tuple
      mass_matrix = ntuple(i -> I, length(f))
    end

    if jac === nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
      if iip
        jac = update_coefficients! # (J, u, p, t)
      else
        jac = (u, p, t) -> update_coefficients!(deepcopy(jac_prototype), u, p, t)
      end
    end

    if recompile
      new{iip,unpack,typeof(f),typeof(mass_matrix),typeof(analytic),typeof(tgrad),
          typeof(jac),typeof(jac_prototype),typeof(Wfact),typeof(Wfact_t),
          typeof(paramjac),typeof(syms),typeof(colorvec)}(
            f,mass_matrix,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
            paramjac,syms,colorvec)
    else
      new{iip,unpack,Any,Any,Any,Any,Any,Any,Any,Any,Any,typeof(syms),typeof(colorvec)}(
        f,mass_matrix,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
        paramjac,syms,colorvec)
    end
  end

  function ODEFunction{iip,recompile,unpack}(f::ODEFunction;
                                             kwargs...) where {iip,recompile,unpack}
    ODEFunction{iip,recompile,unpack}(f.f;
                                      mass_matrix = f.mass_matrix,
                                      analytic = f.analytic,
                                      tgrad = f.tgrad,
                                      jac = f.jac,
                                      jac_prototype = f.jac_prototype,
                                      Wfact = f.Wfact,
                                      Wfact_t = f.Wfact_t,
                                      paramjac = f.paramjac,
                                      syms = f.syms,
                                      colorvec = f.colorvec,
                                      kwargs...)
  end

  ODEFunction{iip,recompile}(f; kwargs...) where {iip,recompile} =
    ODEFunction{iip,RECOMPILE_BY_DEFAULT,UNPACK_BY_DEFAULT}(f; kwargs...)

  ODEFunction{iip}(f; kwargs...) where iip =
    ODEFunction{iip,RECOMPILE_BY_DEFAULT}(f; kwargs...)
end

ODEFunction(f; kwargs...) = ODEFunction{isinplace(f,4)}(f; kwargs...)
ODEFunction(f::ODEFunction; kwargs...) = ODEFunction{isinplace(f)}(f; kwargs...)

"""
$(TYPEDEF)

TODO
"""
struct SplitFunction{iip,unpack,F1,F2,TMM,C,Ta} <: AbstractODEFunction{iip}
  f1::F1
  f2::F2
  mass_matrix::TMM
  cache::C
  analytic::Ta

  function SplitFunction{iip,recompile,unpack}(f1, f2;
                                               mass_matrix = I,
                                               cache = nothing,
                                               analytic = nothing) where {iip,recompile,unpack}
    if recompile
      new{iip,unpack,typeof(f1),typeof(f2),typeof(mass_matrix),typeof(cache),
          typeof(analytic)}(f1,f2,mass_matrix,cache,analytic)
    else
      new{iip,unpack,Any,Any,Any,Any,Any}(f1,f2,mass_matrix,cache,analytic)
    end
  end

  function SplitFunction{iip,recompile,unpack}(f::SplitFunction;
                                               kwargs...) where {iip,recompile,unpack}
    SplitFunction{iip,recompile,unpack}(f.f1, f.f2;
                                        mass_matrix = f.mass_matrix,
                                        cache = f.cache,
                                        analytic = f.analytic,
                                        kwargs...)
  end

  SplitFunction{iip,recompile}(args...; kwargs...) where {iip,recompile} =
    SplitFunction{iip,RECOMPILE_BY_DEFAULT,UNPACK_BY_DEFAULT}(args...; kwargs...)

  SplitFunction{iip}(args...; kwargs...) where iip =
    SplitFunction{iip,RECOMPILE_BY_DEFAULT}(args...; kwargs...)
end

SplitFunction(f1, f2; kwargs...) = SplitFunction{isinplace(f2,4)}(f1, f2; kwargs...)
SplitFunction(f::SplitFunction; kwargs...) = SplitFunction{isinplace(f)}(f; kwargs...)
@add_kwonly SplitFunction(f1, f2, mass_matrix, cache, analytic) =
  SplitFunction(f1, f2; mass_matrix = mass_matrix, cache = cache, analytic = analytic)

"""
$(TYPEDEF)

TODO
"""
struct DynamicalODEFunction{iip,unpack,F1,F2,TMM,Ta} <: AbstractODEFunction{iip}
  f1::F1
  f2::F2
  mass_matrix::TMM
  analytic::Ta

  function DynamicalODEFunction{iip,recompile,unpack}(f1, f2;
                                                      mass_matrix = (I, I),
                                                      analytic = nothing) where {iip,recompile,unpack}
    if recompile
      new{iip,unpack,typeof(f1),typeof(f2),typeof(mass_matrix),typeof(analytic)}(
        f1,f2,mass_matrix,analytic)
    else
      new{iip,unpack,Any,Any,Any,Any}(f1,f2,mass_matrix,analytic)
    end
  end

  function DynamicalODEFunction{iip,recompile,unpack}(f::DynamicalODEFunction;
                                                      kwargs...) where {iip,recompile,unpack}
    DynamicalODEFunction{iip,recompile,unpack}(f.f1, f.f2;
                                               mass_matrix = f.mass_matrix,
                                               analytic = f.analytic,
                                               kwargs...)
  end

  function DynamicalODEFunction{iip,recompile}(args...; kwargs...) where {iip,recompile}
    DynamicalODEFunction{iip,recompile,UNPACK_BY_DEFAULT}(args...; kwargs...)
  end

  function DynamicalODEFunction{iip}(args...; kwargs...) where iip
    DynamicalODEFunction{iip,RECOMPILE_BY_DEFAULT}(args...; kwargs...)
  end
end

DynamicalODEFunction(f1, f2 = nothing; kwargs...) =
  DynamicalODEFunction{isinplace(f1,5)}(f1, f2; kwargs...)
DynamicalODEFunction(f::DynamicalODEFunction; kwargs...) =
  DynamicalODEFunction{isinplace(f)}(f; kwargs...)
@add_kwonly DynamicalODEFunction(f1, f2, mass_matrix, analytic) =
  DynamicalODEFunction(f1, f2; mass_matrix = mass_matrix, analytic = analytic)

"""
$(TYPEDEF)

TODO
"""
abstract type AbstractDDEFunction{iip} <: AbstractDiffEqFunction{iip} end

"""
$(TYPEDEF)

TODO
"""
struct DDEFunction{iip,F,TMM,Ta,Tt,TJ,JP,TW,TWt,TPJ,S,TCV} <: AbstractDDEFunction{iip}
  f::F
  mass_matrix::TMM
  analytic::Ta
  tgrad::Tt
  jac::TJ
  jac_prototype::JP
  Wfact::TW
  Wfact_t::TWt
  paramjac::TPJ
  syms::S
  colorvec::TCV
end

"""
$(TYPEDEF)

TODO
"""
abstract type AbstractDiscreteFunction{iip} <: AbstractDiffEqFunction{iip} end

"""
$(TYPEDEF)

TODO
"""
struct DiscreteFunction{iip,unpack,F,Ta,S} <: AbstractDiscreteFunction{iip}
  f::F
  analytic::Ta
  syms::S

  function DiscreteFunction{iip,recompile,unpack}(f;
                                                  analytic = nothing,
                                                  syms = nothing) where {iip,recompile,unpack}
    if recompile
      new{iip,unpack,typeof(f),typeof(analytic),typeof(syms)}(f,analytic,syms)
    else
      new{iip,unpack,Any,Any,Any}(f,analytic,syms)
    end
  end

  function DiscreteFunction{iip,recompile,unpack}(f::DiscreteFunction;
                                                  kwargs...) where {iip,recompile,unpack}
    DiscreteFunction{iip,recompile,unpack}(f.f; analytic = f.analytic, syms = f.syms, kwargs...)
  end

  DiscreteFunction{iip,recompile}(f; kwargs...) where {iip,recompile} =
    DiscreteFunction{iip,recompile,UNPACK_BY_DEFAULT}(f; kwargs...)

  DiscreteFunction{iip}(f; kwargs...) where iip =
    DiscreteFunction{iip,RECOMPILE_BY_DEFAULT}(f; kwargs...)
end

DiscreteFunction(f; kwargs...) = DiscreteFunction{isinplace(f,4)}(f; kwargs...)
DiscreteFunction(f::DiscreteFunction; kwargs...) = DiscreteFunction{isinplace(f)}(f; kwargs...)

"""
$(TYPEDEF)

TODO
"""
abstract type AbstractSDEFunction{iip} <: AbstractDiffEqFunction{iip} end

"""
$(TYPEDEF)

TODO
"""
struct SDEFunction{iip,F,G,TMM,Ta,Tt,TJ,JP,TW,TWt,TPJ,S,GG,TCV} <: AbstractSDEFunction{iip}
  f::F
  g::G
  mass_matrix::TMM
  analytic::Ta
  tgrad::Tt
  jac::TJ
  jac_prototype::JP
  Wfact::TW
  Wfact_t::TWt
  paramjac::TPJ
  ggprime::GG
  syms::S
  colorvec::TCV
end

"""
$(TYPEDEF)

TODO
"""
struct SplitSDEFunction{iip,F1,F2,G,TMM,C,Ta} <: AbstractSDEFunction{iip}
  f1::F1
  f2::F2
  g::G
  mass_matrix::TMM
  cache::C
  analytic::Ta
end

"""
$(TYPEDEF)

TODO
"""
abstract type AbstractRODEFunction{iip} <: AbstractDiffEqFunction{iip} end

"""
$(TYPEDEF)

TODO
"""
struct RODEFunction{iip,F,TMM,Ta,Tt,TJ,JP,TW,TWt,TPJ,S,TCV} <: AbstractRODEFunction{iip}
  f::F
  mass_matrix::TMM
  analytic::Ta
  tgrad::Tt
  jac::TJ
  jac_prototype::JP
  Wfact::TW
  Wfact_t::TWt
  paramjac::TPJ
  syms::S
  colorvec::TCV
end

"""
$(TYPEDEF)

TODO
"""
abstract type AbstractDAEFunction{iip} <: AbstractDiffEqFunction{iip} end

"""
$(TYPEDEF)

TODO
"""
struct DAEFunction{iip,F,Ta,Tt,TJ,JP,TW,TWt,TPJ,S,TCV} <: AbstractDAEFunction{iip}
  f::F
  analytic::Ta
  tgrad::Tt
  jac::TJ
  jac_prototype::JP
  Wfact::TW
  Wfact_t::TWt
  paramjac::TPJ
  syms::S
  colorvec::TCV
end

######### Functional interface

tgrad(f::ODEFunction{false}, u, p, t) = f.tgrad(u, p, t)
tgrad(f::ODEFunction{false,true}, u, t, integrator::DEIntegrator) =
  tgrad(u, get_params(integrator), t)
tgrad!(f::ODEFunction{true}, dT, u, p, t) = f.tgrad(dT, u, p, t)
tgrad!(f::ODEFunction{true,true}, dT, u, t, integrator::DEIntegrator) =
  tgrad!(dT, u, get_params(integrator), t)

jac(f::ODEFunction{false}, u, p, t) = f.jac(u, p, t)
jac(f::ODEFunction{false,true}, u, t, integrator::DEIntegrator) =
  jac(u, get_params(integrator), t)
jac!(f::ODEFunction{true}, J, u, p, t) = f.jac(J, u, p, t)
jac!(f::ODEFunction{true,true}, J, u, t, integrator::DEIntegrator) =
  jac!(J, u, get_params(integrator), t)

Wfact(f::ODEFunction{false}, u, p, gamma, t) = f.Wfact(u, p, gamma, t)
Wfact(f::ODEFunction{false,true}, u, gamma, t, integrator::DEIntegrator) =
  Wfact(f, u, get_params(integrator), gamma, t)
Wfact!(f::ODEFunction{true}, W, u, p, gamma, t) = f.Wfact(W, u, p, gamma, t)
Wfact!(f::ODEFunction{true,true}, W, u, gamma, t, integrator::DEIntegrator) =
  Wfact!(f, W, u, get_params(integrator), gamma, t)

Wfact_t(f::ODEFunction{false}, u, p, gamma, t) = f.Wfact_t(u, p, gamma, t)
Wfact_t(f::ODEFunction{false,true}, u, gamma, t, integrator::DEIntegrator) =
  Wfact_t(f, u, get_params(integrator), gamma, t)
Wfact_t!(f::ODEFunction{true}, W, u, p, gamma, t) = f.Wfact_t(W, u, p, gamma, t)
Wfact_t!(f::ODEFunction{true,true}, W, u, gamma, t, integrator::DEIntegrator) =
  Wfact_t!(f, W, u, get_params(integrator), gamma, t)

paramjac(f::ODEFunction{false}, u, p, t) = f.paramjac(u, p, t)
paramjac(f::ODEFunction{false,true}, u, t, integrator::DEIntegrator) =
  paramjac(f, u, get_params(integrator), t)
paramjac!(f::ODEFunction{true}, pJ, u, p, t) = f.paramjac(pJ, u, p, t)
paramjac!(f::ODEFunction{true,true}, pJ, u, t, integrator::DEIntegrator) =
  paramjac!(f, pJ, u, get_params(integrator), t)

for F in (:f1, :f2)
  @eval begin
    $F(f::DynamicalODEFunction{false}, u, p, t) = f.$F(u, p, t)
    $F(f::DynamicalODEFunction{false,true}, u, t, integrator::DEIntegrator) =
      $F(f, u, get_params(integrator), t)
    $(Symbol(F, :!))(f!::DynamicalODEFunction{true}, du, u, p, t) = f.$F(du, u, p, t)
    $(Symbol(F, :!))(f!::DynamicalODEFunction{true,true}, du, u, t, integrator::DEIntegrator) =
      $(Symbol(F, :!))(f!, du, u, get_params(integrator), t)
  end
end

######### Function overloads

(f::ODEFunction{false})(u, p, t) = f.f(u, p, t)
(f::ODEFunction{false,true})(u, t, integrator::DEIntegrator) =
  f(u, get_params(integrator), t)

(f!::ODEFunction{true})(du, u, p, t) = f!.f(du, u, p, t)
(f!::ODEFunction{true,true})(du, u, t, integrator::DEIntegrator) =
  f!(du, u, get_params(integrator), t)

(f::DynamicalODEFunction{false})(u, p, t) =
  ArrayPartition(f1(f, u.x[1], u.x[2], p, t), f2(f, u.x[1], u.x[2], p, t))
(f::DynamicalODEFunction{false,true})(u, t, integrator::DEIntegrator) =
  f(u, get_params(integrator), t)
function (f!::DynamicalODEFunction{true})(du, u, p, t)
  f1!(f!, du.x[1], u.x[1], u.x[2], p, t)
  f2!(f!, du.x[2], u.x[1], u.x[2], p, t)
  nothing
end
(f!::DynamicalODEFunction{true,true})(u, t, integrator::DEIntegrator) =
  f!(du, u, get_params(integrator), t)

(f::SplitFunction{false})(u, p, t) = f1(f, u, p, t) + f2(f, u, p, t)
(f::SplitFunction{false,true})(u, t, integrator::DEIntegrator) =
  f(u, get_params(integrator), t)
function (f!::SplitFunction{true})(du, u, p, t)
  f1!(f!, f!.cache, u, p, t)
  f2!(f!, du, u, p, t)
  du .+= f!.cache
  nothing
end
(f!::SplitFunction{true,true})(du, u, t, integrator::DEIntegrator) =
  f!(du, u, get_params(integrator), t)

(f::DiscreteFunction{false})(u, p, t) = f.f(u, p, t)
(f::DiscreteFunction{false,true})(u, t, integrator::DEIntegrator) =
  f(u, get_params(integrator), t)
(f!::DiscreteFunction{true})(du, u, p, t) = f!.f(du, u, p, t)
(f!::DiscreteFunction{true,true})(du, u, t, integrator::DEIntegrator) =
  f!(du, u, get_params(integrator), t)

######### Backwards Compatibility Overloads

(f::DAEFunction)(args...) = f.f(args...)
(f::DAEFunction)(::Type{Val{:analytic}},args...) = f.analytic(args...)
(f::DAEFunction)(::Type{Val{:tgrad}},args...) = f.tgrad(args...)
(f::DAEFunction)(::Type{Val{:jac}},args...) = f.jac(args...)
(f::DAEFunction)(::Type{Val{:Wfact}},args...) = f.Wfact(args...)
(f::DAEFunction)(::Type{Val{:Wfact_t}},args...) = f.Wfact_t(args...)
(f::DAEFunction)(::Type{Val{:paramjac}},args...) = f.paramjac(args...)

(f::DDEFunction)(args...) = f.f(args...)
(f::DDEFunction)(::Type{Val{:analytic}},args...) = f.analytic(args...)

(f::SDEFunction)(args...) = f.f(args...)
(f::SDEFunction)(::Type{Val{:analytic}},args...) = f.analytic(args...)
(f::SDEFunction)(::Type{Val{:tgrad}},args...) = f.tgrad(args...)
(f::SDEFunction)(::Type{Val{:jac}},args...) = f.jac(args...)
(f::SDEFunction)(::Type{Val{:Wfact}},args...) = f.Wfact(args...)
(f::SDEFunction)(::Type{Val{:Wfact_t}},args...) = f.Wfact_t(args...)
(f::SDEFunction)(::Type{Val{:paramjac}},args...) = f.paramjac(args...)

(f::SplitSDEFunction)(u,p,t) = f.f1(u,p,t) + f.f2(u,p,t)
(f::SplitSDEFunction)(::Type{Val{:analytic}},args...) = f.analytic(args...)
function (f::SplitSDEFunction)(du,u,p,t)
  f.f1(f.cache,u,p,t)
  f.f2(du,u,p,t)
  du .+= f.cache
end

(f::RODEFunction)(args...) = f.f(args...)

######### Basic Constructor

function SDEFunction{iip,true}(f,g;
                 mass_matrix=I,
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 jac_prototype=nothing,
                 Wfact=nothing,
                 Wfact_t=nothing,
                 paramjac = nothing,
                 ggprime = nothing,
                 syms = nothing,
                 colorvec = nothing) where iip
                 if jac == nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                  if iip
                    jac = update_coefficients! #(J,u,p,t)
                  else
                    jac = (u,p,t) -> update_coefficients!(deepcopy(jac_prototype),u,p,t)
                  end
                 end
                 SDEFunction{iip,typeof(f),typeof(g),
                 typeof(mass_matrix),typeof(analytic),typeof(tgrad),
                 typeof(jac),typeof(jac_prototype),typeof(Wfact),typeof(Wfact_t),
                 typeof(paramjac),typeof(syms),
                 typeof(ggprime),typeof(colorvec)}(
                 f,g,mass_matrix,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
                 paramjac,ggprime,syms,colorvec)
end
function SDEFunction{iip,false}(f,g;
                 mass_matrix=I,
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 jac_prototype=nothing,
                 Wfact=nothing,
                 Wfact_t=nothing,
                 paramjac = nothing,
                 ggprime = nothing,
                 syms = nothing,
                 colorvec = nothing) where iip
                 if jac == nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                  if iip
                    jac = update_coefficients! #(J,u,p,t)
                  else
                    jac = (u,p,t) -> update_coefficients!(deepcopy(jac_prototype),u,p,t)
                  end
                 end
                 SDEFunction{iip,Any,Any,Any,Any,Any,
                 Any,Any,Any,Any,
                 Any,typeof(syms),Any,typeof(colorvec)}(
                 f,g,mass_matrix,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
                 paramjac,ggprime,syms,colorvec)
end
SDEFunction{iip}(f,g; kwargs...) where iip = SDEFunction{iip,RECOMPILE_BY_DEFAULT}(f,g; kwargs...)
SDEFunction{iip}(f::SDEFunction,g; kwargs...) where iip = f
SDEFunction(f,g; kwargs...) = SDEFunction{isinplace(f, 4),RECOMPILE_BY_DEFAULT}(f,g; kwargs...)
SDEFunction(f::SDEFunction; kwargs...) = f

SplitSDEFunction{iip,true}(f1,f2,g; mass_matrix=I,
                           _func_cache=nothing,analytic=nothing) where iip =
SplitSDEFunction{iip,typeof(f1),typeof(f2),typeof(g),
              typeof(mass_matrix),typeof(_func_cache),
              typeof(analytic)}(f1,f2,g,mass_matrix,_func_cache,analytic)
SplitSDEFunction{iip,false}(f1,f2,g; mass_matrix=I,
                            _func_cache=nothing,analytic=nothing) where iip =
SplitSDEFunction{iip,Any,Any,Any,Any,Any}(f1,f2,g,mass_matrix,_func_cache,analytic)
SplitSDEFunction(f1,f2,g; kwargs...) = SplitSDEFunction{isinplace(f2, 4)}(f1, f2, g; kwargs...)
SplitSDEFunction{iip}(f1,f2, g; kwargs...) where iip =
SplitSDEFunction{iip,RECOMPILE_BY_DEFAULT}(SDEFunction(f1,g), SDEFunction{iip}(f2,g), g; kwargs...)
SplitSDEFunction(f::SplitSDEFunction; kwargs...) = f

function RODEFunction{iip,true}(f;
                 mass_matrix=I,
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 jac_prototype=nothing,
                 Wfact=nothing,
                 Wfact_t=nothing,
                 paramjac = nothing,
                 syms = nothing,
                 colorvec = nothing) where iip
                 if jac == nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                  if iip
                    jac = update_coefficients! #(J,u,p,t)
                  else
                    jac = (u,p,t) -> update_coefficients!(deepcopy(jac_prototype),u,p,t)
                  end
                 end
                 RODEFunction{iip,typeof(f),typeof(mass_matrix),
                 typeof(analytic),typeof(tgrad),
                 typeof(jac),typeof(jac_prototype),typeof(Wfact),typeof(Wfact_t),
                 typeof(paramjac),typeof(syms),typeof(colorvec)}(
                 f,mass_matrix,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
                 paramjac,syms,colorvec)
end
function RODEFunction{iip,false}(f;
                 mass_matrix=I,
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 Wfact=nothing,
                 Wfact_t=nothing,
                 paramjac = nothing,
                 syms = nothing,
                 colorvec = nothing) where iip
                 if jac == nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                  if iip
                    jac = update_coefficients! #(J,u,p,t)
                  else
                    jac = (u,p,t) -> update_coefficients!(deepcopy(jac_prototype),u,p,t)
                  end
                 end
                 RODEFunction{iip,Any,Any,Any,Any,
                 Any,Any,Any,Any,
                 Any,typeof(syms),typeof(colorvec)}(
                 f,mass_matrix,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
                 paramjac,syms,colorvec)
end
RODEFunction{iip}(f; kwargs...) where iip = RODEFunction{iip,RECOMPILE_BY_DEFAULT}(f; kwargs...)
RODEFunction{iip}(f::RODEFunction; kwargs...) where iip = f
RODEFunction(f; kwargs...) = RODEFunction{isinplace(f, 5),RECOMPILE_BY_DEFAULT}(f; kwargs...)
RODEFunction(f::RODEFunction; kwargs...) = f

function DAEFunction{iip,true}(f;
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 jac_prototype=nothing,
                 Wfact=nothing,
                 Wfact_t=nothing,
                 paramjac = nothing,
                 syms = nothing,
                 colorvec = nothing) where iip
                 if jac == nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                  if iip
                    jac = update_coefficients! #(J,u,p,t)
                  else
                    jac = (u,p,t) -> update_coefficients!(deepcopy(jac_prototype),u,p,t)
                  end
                 end
                 DAEFunction{iip,typeof(f),typeof(analytic),typeof(tgrad),
                 typeof(jac),typeof(jac_prototype),typeof(Wfact),typeof(Wfact_t),
                 typeof(paramjac),typeof(syms),typeof(colorvec)}(
                 f,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
                 paramjac,syms,colorvec)
end
function DAEFunction{iip,false}(f;
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 jac_prototype=nothing,
                 Wfact=nothing,
                 Wfact_t=nothing,
                 paramjac = nothing,
                 syms = nothing,
                 colorvec = nothing) where iip
                 if jac == nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                  if iip
                    jac = update_coefficients! #(J,u,p,t)
                  else
                    jac = (u,p,t) -> update_coefficients!(deepcopy(jac_prototype),u,p,t)
                  end
                 end
                 DAEFunction{iip,Any,Any,Any,
                 Any,Any,Any,Any,
                 Any,typeof(syms),typeof(colorvec)}(
                 f,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
                 paramjac,syms,colorvec)
end
DAEFunction{iip}(f; kwargs...) where iip = DAEFunction{iip,RECOMPILE_BY_DEFAULT}(f; kwargs...)
DAEFunction{iip}(f::DAEFunction; kwargs...) where iip = f
DAEFunction(f; kwargs...) = DAEFunction{isinplace(f, 5),RECOMPILE_BY_DEFAULT}(f; kwargs...)
DAEFunction(f::DAEFunction; kwargs...) = f

function DDEFunction{iip,true}(f;
                 mass_matrix=I,
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 jac_prototype=nothing,
                 Wfact=nothing,
                 Wfact_t=nothing,
                 paramjac = nothing,
                 syms = nothing,
                 colorvec = nothing) where iip
                 if jac == nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                  if iip
                    jac = update_coefficients! #(J,u,p,t)
                  else
                    jac = (u,p,t) -> update_coefficients!(deepcopy(jac_prototype),u,p,t)
                  end
                 end
                 DDEFunction{iip,typeof(f),typeof(mass_matrix),typeof(analytic),typeof(tgrad),
                 typeof(jac),typeof(jac_prototype),typeof(Wfact),typeof(Wfact_t),
                 typeof(paramjac),typeof(syms),typeof(colorvec)}(
                 f,mass_matrix,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
                 paramjac,syms,colorvec)
end
function DDEFunction{iip,false}(f;
                 mass_matrix=I,
                 analytic=nothing,
                 tgrad=nothing,
                 jac=nothing,
                 jac_prototype=nothing,
                 Wfact=nothing,
                 Wfact_t=nothing,
                 paramjac = nothing,
                 syms = nothing,
                 colorvec = nothing) where iip
                 if jac == nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                  if iip
                    jac = update_coefficients! #(J,u,p,t)
                  else
                    jac = (u,p,t) -> update_coefficients!(deepcopy(jac_prototype),u,p,t)
                  end
                 end
                 DDEFunction{iip,Any,Any,Any,Any,
                 Any,Any,Any,Any,
                 Any,typeof(syms),typeof(colorvec)}(
                 f,mass_matrix,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
                 paramjac,syms,colorvec)
end
DDEFunction{iip}(f; kwargs...) where iip = DDEFunction{iip,RECOMPILE_BY_DEFAULT}(f; kwargs...)
DDEFunction{iip}(f::DDEFunction; kwargs...) where iip = f
DDEFunction(f; kwargs...) = DDEFunction{isinplace(f, 5),RECOMPILE_BY_DEFAULT}(f; kwargs...)
DDEFunction(f::DDEFunction; kwargs...) = f

########## Existance Functions

# compatibility
has_invW(f::AbstractDiffEqFunction) = false
has_analytic(f::AbstractDiffEqFunction) = f.analytic != nothing
has_jac(f::AbstractDiffEqFunction) = f.jac != nothing
has_tgrad(f::AbstractDiffEqFunction) = f.tgrad != nothing
has_Wfact(f::AbstractDiffEqFunction) = f.Wfact != nothing
has_Wfact_t(f::AbstractDiffEqFunction) = f.Wfact_t != nothing
has_paramjac(f::AbstractDiffEqFunction) = f.paramjac != nothing
has_syms(f::AbstractDiffEqFunction) = :syms ∈ fieldnames(typeof(f)) && f.syms != nothing
has_colorvec(f::AbstractDiffEqFunction) = f.colorvec != nothing

# TODO: find an appropriate way to check `has_*`
has_jac(f::Union{SplitFunction,SplitSDEFunction}) = has_jac(f.f1)
has_tgrad(f::Union{SplitFunction,SplitSDEFunction}) = has_tgrad(f.f1)
has_Wfact(f::Union{SplitFunction,SplitSDEFunction}) = has_Wfact(f.f1)
has_Wfact_t(f::Union{SplitFunction,SplitSDEFunction}) = has_Wfact_t(f.f1)
has_paramjac(f::Union{SplitFunction,SplitSDEFunction}) = has_paramjac(f.f1)
has_colorvec(f::Union{SplitFunction,SplitSDEFunction}) = has_colorvec(f.f1)

has_jac(f::DynamicalODEFunction) = has_jac(f.f1)
has_tgrad(f::DynamicalODEFunction) = has_tgrad(f.f1)
has_Wfact(f::DynamicalODEFunction) = has_Wfact(f.f1)
has_Wfact_t(f::DynamicalODEFunction) = has_Wfact_t(f.f1)
has_paramjac(f::DynamicalODEFunction) = has_paramjac(f.f1)
has_colorvec(f::DynamicalODEFunction) = has_colorvec(f.f1)

######### Additional traits

islinear(f) = false # fallback
islinear(::AbstractDiffEqFunction) = false
islinear(f::ODEFunction) = islinear(f.f)
islinear(f::SplitFunction) = islinear(f.f1)

######### Compatibility Constructor from Tratis

function Base.convert(::Type{ODEFunction}, f)
  if __has_analytic(f)
    analytic = f.analytic
  else
    analytic = nothing
  end
  if __has_jac(f)
    jac = f.jac
  else
    jac = nothing
  end
  if __has_tgrad(f)
    tgrad = f.tgrad
  else
    tgrad = nothing
  end
  if __has_Wfact(f)
    Wfact = f.Wfact
  else
    Wfact = nothing
  end
  if __has_Wfact_t(f)
    Wfact_t = f.Wfact_t
  else
    Wfact_t = nothing
  end
  if __has_paramjac(f)
    paramjac = f.paramjac
  else
    paramjac = nothing
  end
  if __has_syms(f)
    syms = f.syms
  else
    syms = nothing
  end
  if __has_colorvec(f)
    colorvec = f.colorvec
  else
    colorvec = nothing
  end
  ODEFunction(f;analytic=analytic,tgrad=tgrad,jac=jac,Wfact=Wfact,
              Wfact_t=Wfact_t,paramjac=paramjac,syms=syms,colorvec=colorvec)
end
function Base.convert(::Type{ODEFunction{iip}},f) where iip
  if __has_analytic(f)
    analytic = f.analytic
  else
    analytic = nothing
  end
  if __has_jac(f)
    jac = f.jac
  else
    jac = nothing
  end
  if __has_tgrad(f)
    tgrad = f.tgrad
  else
    tgrad = nothing
  end
  if __has_Wfact(f)
    Wfact = f.Wfact
  else
    Wfact = nothing
  end
  if __has_Wfact_t(f)
    Wfact_t = f.Wfact_t
  else
    Wfact_t = nothing
  end
  if __has_paramjac(f)
    paramjac = f.paramjac
  else
    paramjac = nothing
  end
  if __has_syms(f)
    syms = f.syms
  else
    syms = nothing
  end
  if __has_colorvec(f)
    colorvec = f.colorvec
  else
    colorvec = nothing
  end
  ODEFunction{iip,RECOMPILE_BY_DEFAULT}(f;analytic=analytic,tgrad=tgrad,jac=jac,Wfact=Wfact,
              Wfact_t=Wfact_t,paramjac=paramjac,syms=syms,colorvec=colorvec)
end

function Base.convert(::Type{DiscreteFunction},f)
  if __has_analytic(f)
    analytic = f.analytic
  else
    analytic = nothing
  end
  if __has_syms(f)
    syms = f.syms
  else
    syms = nothing
  end
  DiscreteFunction(f;analytic=analytic,syms=syms)
end
function Base.convert(::Type{DiscreteFunction{iip}},f) where iip
  if __has_analytic(f)
    analytic = f.analytic
  else
    analytic = nothing
  end
  if __has_syms(f)
    syms = f.syms
  else
    syms = nothing
  end
  DiscreteFunction{iip,RECOMPILE_BY_DEFAULT}(f;analytic=analytic,syms=syms)
end

function Base.convert(::Type{DAEFunction},f)
  if __has_analytic(f)
    analytic = f.analytic
  else
    analytic = nothing
  end
  if __has_jac(f)
    jac = f.jac
  else
    jac = nothing
  end
  if __has_tgrad(f)
    tgrad = f.tgrad
  else
    tgrad = nothing
  end
  if __has_Wfact(f)
    Wfact = f.Wfact
  else
    Wfact = nothing
  end
  if __has_Wfact_t(f)
    Wfact_t = f.Wfact_t
  else
    Wfact_t = nothing
  end
  if __has_paramjac(f)
    paramjac = f.paramjac
  else
    paramjac = nothing
  end
  if __has_syms(f)
    syms = f.syms
  else
    syms = nothing
  end
  if __has_colorvec(f)
    colorvec = f.colorvec
  else
    colorvec = nothing
  end
  DAEFunction(f;analytic=analytic,tgrad=tgrad,jac=jac,Wfact=Wfact,
              Wfact_t=Wfact_t,paramjac=paramjac,syms=syms,colorvec=colorvec)
end
function Base.convert(::Type{DAEFunction{iip}},f) where iip
  if __has_analytic(f)
    analytic = f.analytic
  else
    analytic = nothing
  end
  if __has_jac(f)
    jac = f.jac
  else
    jac = nothing
  end
  if __has_tgrad(f)
    tgrad = f.tgrad
  else
    tgrad = nothing
  end
  if __has_Wfact(f)
    Wfact = f.Wfact
  else
    Wfact = nothing
  end
  if __has_Wfact_t(f)
    Wfact_t = f.Wfact_t
  else
    Wfact_t = nothing
  end
  if __has_paramjac(f)
    paramjac = f.paramjac
  else
    paramjac = nothing
  end
  if __has_syms(f)
    syms = f.syms
  else
    syms = nothing
  end
  if __has_colorvec(f)
    colorvec = f.colorvec
  else
    colorvec = nothing
  end
  DAEFunction{iip,RECOMPILE_BY_DEFAULT}(f;analytic=analytic,tgrad=tgrad,jac=jac,Wfact=Wfact,
              Wfact_t=Wfact_t,paramjac=paramjac,syms=syms,colorvec=colorvec)
end

function Base.convert(::Type{DDEFunction},f)
  if __has_analytic(f)
    analytic = f.analytic
  else
    analytic = nothing
  end
  if __has_syms(f)
    syms = f.syms
  else
    syms = nothing
  end
  if __has_colorvec(f)
    colorvec = f.colorvec
  else
    colorvec = nothing
  end
  DDEFunction(f;analytic=analytic,syms=syms,colorvec=colorvec)
end
function Base.convert(::Type{DDEFunction{iip}},f) where iip
  if __has_analytic(f)
    analytic = f.analytic
  else
    analytic = nothing
  end
  if __has_syms(f)
    syms = f.syms
  else
    syms = nothing
  end
  if __has_colorvec(f)
    colorvec = f.colorvec
  else
    colorvec = nothing
  end
  DDEFunction{iip,RECOMPILE_BY_DEFAULT}(f;analytic=analytic,syms=syms,colorvec=colorvec)
end

function Base.convert(::Type{SDEFunction},f,g)
  if __has_analytic(f)
    analytic = f.analytic
  else
    analytic = nothing
  end
  if __has_jac(f)
    jac = f.jac
  else
    jac = nothing
  end
  if __has_tgrad(f)
    tgrad = f.tgrad
  else
    tgrad = nothing
  end
  if __has_Wfact(f)
    Wfact = f.Wfact
  else
    Wfact = nothing
  end
  if __has_Wfact_t(f)
    Wfact_t = f.Wfact_t
  else
    Wfact_t = nothing
  end
  if __has_paramjac(f)
    paramjac = f.paramjac
  else
    paramjac = nothing
  end
  if __has_syms(f)
    syms = f.syms
  else
    syms = nothing
  end
  if __has_colorvec(f)
    colorvec = f.colorvec
  else
    colorvec = nothing
  end
  SDEFunction(f,g;analytic=analytic,tgrad=tgrad,jac=jac,Wfact=Wfact,
              Wfact_t=Wfact_t,paramjac=paramjac,syms=syms,colorvec=colorvec)
end
function Base.convert(::Type{SDEFunction{iip}},f,g) where iip
  if __has_analytic(f)
    analytic = f.analytic
  else
    analytic = nothing
  end
  if __has_jac(f)
    jac = f.jac
  else
    jac = nothing
  end
  if __has_tgrad(f)
    tgrad = f.tgrad
  else
    tgrad = nothing
  end
  if __has_Wfact(f)
    Wfact = f.Wfact
  else
    Wfact = nothing
  end
  if __has_Wfact_t(f)
    Wfact_t = f.Wfact_t
  else
    Wfact_t = nothing
  end
  if __has_paramjac(f)
    paramjac = f.paramjac
  else
    paramjac = nothing
  end
  if __has_syms(f)
    syms = f.syms
  else
    syms = nothing
  end
  if __has_colorvec(f)
    colorvec = f.colorvec
  else
    colorvec = nothing
  end
  SDEFunction{iip,RECOMPILE_BY_DEFAULT}(f,g;analytic=analytic,
              tgrad=tgrad,jac=jac,Wfact=Wfact,
              Wfact_t=Wfact_t,paramjac=paramjac,syms=syms,colorvec=colorvec)
end

function Base.convert(::Type{RODEFunction},f)
  if __has_analytic(f)
    analytic = f.analytic
  else
    analytic = nothing
  end
  if __has_jac(f)
    jac = f.jac
  else
    jac = nothing
  end
  if __has_tgrad(f)
    tgrad = f.tgrad
  else
    tgrad = nothing
  end
  if __has_Wfact(f)
    Wfact = f.Wfact
  else
    Wfact = nothing
  end
  if __has_Wfact_t(f)
    Wfact_t = f.Wfact_t
  else
    Wfact_t = nothing
  end
  if __has_paramjac(f)
    paramjac = f.paramjac
  else
    paramjac = nothing
  end
  if __has_syms(f)
    syms = f.syms
  else
    syms = nothing
  end
  if __has_colorvec(f)
    colorvec = f.colorvec
  else
    colorvec = nothing
  end
  RODEFunction(f;analytic=analytic,tgrad=tgrad,jac=jac,Wfact=Wfact,
              Wfact_t=Wfact_t,paramjac=paramjac,syms=syms,colorvec=colorvec)
end
function Base.convert(::Type{RODEFunction{iip}},f) where iip
  if __has_analytic(f)
    analytic = f.analytic
  else
    analytic = nothing
  end
  if __has_jac(f)
    jac = f.jac
  else
    jac = nothing
  end
  if __has_tgrad(f)
    tgrad = f.tgrad
  else
    tgrad = nothing
  end
  if __has_Wfact(f)
    Wfact = f.Wfact
  else
    Wfact = nothing
  end
  if __has_Wfact_t(f)
    Wfact_t = f.Wfact_t
  else
    Wfact_t = nothing
  end
  if __has_paramjac(f)
    paramjac = f.paramjac
  else
    paramjac = nothing
  end
  if __has_syms(f)
    syms = f.syms
  else
    syms = nothing
  end
  if __has_colorvec(f)
    colorvec = f.colorvec
  else
    colorvec = nothing
  end
  RODEFunction{iip,RECOMPILE_BY_DEFAULT}(f;analytic=analytic,
              tgrad=tgrad,jac=jac,Wfact=Wfact,
              Wfact_t=Wfact_t,paramjac=paramjac,syms=syms,colorvec=colorvec)
end
