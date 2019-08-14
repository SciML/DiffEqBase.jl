const RECOMPILE_BY_DEFAULT = true

"""
$(TYPEDEF)

TODO
"""
abstract type AbstractODEFunction{iip} <: AbstractDiffEqFunction{iip} end

"""
$(TYPEDEF)

TODO
"""
struct ODEFunction{iip,F,TMM,Ta,Tt,TJ,JP,TW,TWt,TPJ,S,TCV} <: AbstractODEFunction{iip}
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
struct SplitFunction{iip,F1,F2,TMM,C,Ta,Tt,TJ,JP,TW,TWt,TPJ,S,TCV} <: AbstractODEFunction{iip}
  f1::F1
  f2::F2
  mass_matrix::TMM
  cache::C
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
struct DynamicalODEFunction{iip,F1,F2,TMM,Ta} <: AbstractODEFunction{iip}
  f1::F1
  f2::F2
  mass_matrix::TMM
  analytic::Ta
end

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
struct DiscreteFunction{iip,F,Ta,S} <: AbstractDiscreteFunction{iip}
  f::F
  analytic::Ta
  syms::S
end

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
struct SplitSDEFunction{iip,F1,F2,G,TMM,C,Ta,Tt,TJ,JP,TW,TWt,TPJ,S,TCV} <: AbstractSDEFunction{iip}
  f1::F1
  f2::F2
  g::G
  mass_matrix::TMM
  cache::C
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

######### Backwards Compatibility Overloads

(f::ODEFunction)(args...) = f.f(args...)
(f::ODEFunction)(::Type{Val{:analytic}},args...) = f.analytic(args...)
(f::ODEFunction)(::Type{Val{:tgrad}},args...) = f.tgrad(args...)
(f::ODEFunction)(::Type{Val{:jac}},args...) = f.jac(args...)
(f::ODEFunction)(::Type{Val{:Wfact}},args...) = f.Wfact(args...)
(f::ODEFunction)(::Type{Val{:Wfact_t}},args...) = f.Wfact_t(args...)
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

function ODEFunction{iip,true}(f;
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

                 if mass_matrix == I && typeof(f) <: Tuple
                  mass_matrix = ((I for i in 1:length(f))...,)
                 end

                 if jac === nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                  if iip
                    jac = update_coefficients! #(J,u,p,t)
                  else
                    jac = (u,p,t) -> update_coefficients!(deepcopy(jac_prototype),u,p,t)
                  end
                 end

                 if jac_prototype !== nothing && colorvec === nothing && ArrayInterface.fast_matrix_colors(jac_prototype)
                   _colorvec = ArrayInterface.matrix_colors(jac_prototype)
                 else
                   _colorvec = colorvec
                 end

                 ODEFunction{iip,typeof(f),typeof(mass_matrix),typeof(analytic),typeof(tgrad),
                 typeof(jac),typeof(jac_prototype),typeof(Wfact),typeof(Wfact_t),
                 typeof(paramjac),typeof(syms),typeof(_colorvec)}(
                 f,mass_matrix,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
                 paramjac,syms,_colorvec)
end
function ODEFunction{iip,false}(f;
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

                 if jac === nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                  if iip
                    jac = update_coefficients! #(J,u,p,t)
                  else
                    jac = (u,p,t) -> update_coefficients!(deepcopy(jac_prototype),u,p,t)
                  end
                 end

                 if jac_prototype !== nothing && colorvec === nothing && ArrayInterface.fast_matrix_colors(jac_prototype)
                   _colorvec = ArrayInterface.matrix_colors(jac_prototype)
                 else
                   _colorvec = colorvec
                 end

                 ODEFunction{iip,Any,Any,Any,Any,
                 Any,Any,Any,Any,
                 Any,typeof(syms),typeof(_colorvec)}(
                 f,mass_matrix,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
                 paramjac,syms,_colorvec)
end
ODEFunction{iip}(f; kwargs...) where iip = ODEFunction{iip,RECOMPILE_BY_DEFAULT}(f; kwargs...)
ODEFunction{iip}(f::ODEFunction; kwargs...) where iip = f
ODEFunction(f; kwargs...) = ODEFunction{isinplace(f, 4),RECOMPILE_BY_DEFAULT}(f; kwargs...)
ODEFunction(f::ODEFunction; kwargs...) = f

@add_kwonly function SplitFunction(f1,f2,mass_matrix,cache,analytic,tgrad,jac,
                                   jac_prototype,Wfact,Wfact_t,paramjac,
                                   syms,colorvec)
  f1 = typeof(f1) <: AbstractDiffEqOperator ? f1 : ODEFunction(f1)
  f2 = ODEFunction(f2)
  SplitFunction{isinplace(f2),typeof(f1),typeof(f2),typeof(mass_matrix),
              typeof(cache),typeof(analytic),typeof(tgrad),typeof(jac),
              typeof(jac_prototype),
              typeof(Wfact),typeof(Wfact_t),typeof(paramjac),typeof(syms),
              typeof(colorvec)}(f1,f2,mass_matrix,cache,analytic,tgrad,jac,
              jac_prototype,Wfact,Wfact_t,paramjac,syms,colorvec)
end
function SplitFunction{iip,true}(f1,f2;
                                 mass_matrix=I,_func_cache=nothing,
                                 analytic=nothing,
                                 tgrad = nothing,
                                 jac = nothing,
                                 jac_prototype = nothing,
                                 Wfact = nothing,
                                 Wfact_t = nothing,
                                 paramjac = nothing,
                                 syms = nothing,
                                 colorvec = nothing) where iip
  SplitFunction{iip,typeof(f1),typeof(f2),typeof(mass_matrix),
                typeof(_func_cache),typeof(analytic),
                typeof(tgrad),typeof(jac),typeof(jac_prototype),
                typeof(Wfact),typeof(Wfact_t),typeof(paramjac),typeof(syms),
                typeof(colorvec)}(
                f1,f2,mass_matrix,_func_cache,analytic,tgrad,jac,jac_prototype,
                Wfact,Wfact_t,paramjac,syms,colorvec)
end
function SplitFunction{iip,false}(f1,f2; mass_matrix=I,
                                  _func_cache=nothing,analytic=nothing,
                                  tgrad = nothing,
                                  jac = nothing,
                                  jac_prototype = nothing,
                                  Wfact = nothing,
                                  Wfact_t = nothing,
                                  paramjac = nothing,
                                  syms = nothing,
                                  colorvec = nothing) where iip
  SplitFunction{iip,Any,Any,Any,Any,Any,Any,Any,
                Any,Any,Any,Any,Any}(
                f1,f2,mass_matrix,_func_cache,analytic,tgrad,jac,jac_prototype,
                Wfact,Wfact_t,paramjac,syms,colorvec)
end
SplitFunction(f1,f2; kwargs...) = SplitFunction{isinplace(f2, 4)}(f1, f2; kwargs...)
SplitFunction{iip}(f1,f2; kwargs...) where iip =
SplitFunction{iip,RECOMPILE_BY_DEFAULT}(ODEFunction(f1),ODEFunction{iip}(f2); kwargs...)
SplitFunction(f::SplitFunction; kwargs...) = f

@add_kwonly function DynamicalODEFunction{iip}(f1,f2,mass_matrix,analytic) where iip
  f1 = ODEFunction(f1)
  f2 != nothing && (f2 = ODEFunction(f2))
  DynamicalODEFunction{iip,typeof(f1),typeof(f2),typeof(mass_matrix),typeof(analytic)}(f1,f2,mass_matrix,analytic)
end
DynamicalODEFunction{iip,true}(f1,f2;mass_matrix=(I,I),analytic=nothing) where iip =
DynamicalODEFunction{iip,typeof(f1),typeof(f2),typeof(mass_matrix),typeof(analytic)}(f1,f2,mass_matrix,analytic)
DynamicalODEFunction{iip,false}(f1,f2;mass_matrix=(I,I),analytic=nothing) where iip =
DynamicalODEFunction{iip,Any,Any,Any,Any}(f1,f2,mass_matrix,analytic)
DynamicalODEFunction(f1,f2=nothing; kwargs...) = DynamicalODEFunction{isinplace(f1, 5)}(f1, f2; kwargs...)
DynamicalODEFunction{iip}(f1,f2; kwargs...) where iip =
DynamicalODEFunction{iip,RECOMPILE_BY_DEFAULT}(ODEFunction{iip}(f1), ODEFunction{iip}(f2); kwargs...)
DynamicalODEFunction(f::DynamicalODEFunction; kwargs...) = f

function DiscreteFunction{iip,true}(f;
                 analytic=nothing, syms=nothing) where iip
                 DiscreteFunction{iip,typeof(f),typeof(analytic),typeof(syms)}(
                 f,analytic,syms)
end
function DiscreteFunction{iip,false}(f;
                 analytic=nothing, syms=nothing) where iip
                 DiscreteFunction{iip,Any,Any,Any}(
                 f,analytic,syms)
end
DiscreteFunction{iip}(f; kwargs...) where iip = DiscreteFunction{iip,RECOMPILE_BY_DEFAULT}(f; kwargs...)
DiscreteFunction{iip}(f::DiscreteFunction; kwargs...) where iip = f
DiscreteFunction(f; kwargs...) = DiscreteFunction{isinplace(f, 4),RECOMPILE_BY_DEFAULT}(f; kwargs...)
DiscreteFunction(f::DiscreteFunction; kwargs...) = f

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

                 if jac === nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                  if iip
                    jac = update_coefficients! #(J,u,p,t)
                  else
                    jac = (u,p,t) -> update_coefficients!(deepcopy(jac_prototype),u,p,t)
                  end
                 end

                 if jac_prototype !== nothing && colorvec === nothing && ArrayInterface.fast_matrix_colors(jac_prototype)
                   _colorvec = ArrayInterface.matrix_colors(jac_prototype)
                 else
                   _colorvec = colorvec
                 end

                 SDEFunction{iip,typeof(f),typeof(g),
                 typeof(mass_matrix),typeof(analytic),typeof(tgrad),
                 typeof(jac),typeof(jac_prototype),typeof(Wfact),typeof(Wfact_t),
                 typeof(paramjac),typeof(syms),
                 typeof(ggprime),typeof(_colorvec)}(
                 f,g,mass_matrix,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
                 paramjac,ggprime,syms,_colorvec)
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

                 if jac === nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                  if iip
                    jac = update_coefficients! #(J,u,p,t)
                  else
                    jac = (u,p,t) -> update_coefficients!(deepcopy(jac_prototype),u,p,t)
                  end
                 end

                 if jac_prototype !== nothing && colorvec === nothing && ArrayInterface.fast_matrix_colors(jac_prototype)
                   _colorvec = ArrayInterface.matrix_colors(jac_prototype)
                 else
                   _colorvec = colorvec
                 end

                 SDEFunction{iip,Any,Any,Any,Any,Any,
                 Any,Any,Any,Any,
                 Any,typeof(syms),Any,typeof(_colorvec)}(
                 f,g,mass_matrix,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
                 paramjac,ggprime,syms,_colorvec)
end
SDEFunction{iip}(f,g; kwargs...) where iip = SDEFunction{iip,RECOMPILE_BY_DEFAULT}(f,g; kwargs...)
SDEFunction{iip}(f::SDEFunction,g; kwargs...) where iip = f
SDEFunction(f,g; kwargs...) = SDEFunction{isinplace(f, 4),RECOMPILE_BY_DEFAULT}(f,g; kwargs...)
SDEFunction(f::SDEFunction; kwargs...) = f

@add_kwonly function SplitSDEFunction(f1,f2,g,mass_matrix,cache,analytic,tgrad,jac,
                                   jac_prototype,Wfact,Wfact_t,paramjac,
                                   syms,colorvec)
  f1 = typeof(f1) <: AbstractDiffEqOperator ? f1 : SDEFunction(f1)
  f2 = SDEFunction(f2)
  SplitFunction{isinplace(f2),typeof(f1),typeof(f2),typeof(g),typeof(mass_matrix),
              typeof(cache),typeof(analytic),typeof(tgrad),typeof(jac),
              typeof(Wfact),typeof(Wfact_t),typeof(paramjac),typeof(syms),
              typeof(colorvec)}(f1,f2,mass_matrix,cache,analytic,tgrad,jac,
              jac_prototype,Wfact,Wfact_t,paramjac,syms,colorvec)
end

function SplitSDEFunction{iip,true}(f1,f2,g; mass_matrix=I,
                           _func_cache=nothing,analytic=nothing,
                           tgrad = nothing,
                           jac = nothing,
                           jac_prototype = nothing,
                           Wfact = nothing,
                           Wfact_t = nothing,
                           paramjac = nothing,
                           syms = nothing,
                           colorvec = nothing) where iip
  SplitSDEFunction{iip,typeof(f1),typeof(f2),typeof(g),
              typeof(mass_matrix),typeof(_func_cache),
              typeof(analytic),
              typeof(tgrad),typeof(jac),typeof(jac_prototype),
              typeof(Wfact),typeof(Wfact_t),typeof(paramjac),typeof(syms),
              typeof(colorvec)}(f1,f2,g,mass_matrix,_func_cache,analytic,
              tgrad,jac,jac_prototype,
              Wfact,Wfact_t,paramjac,syms,colorvec)
end
function SplitSDEFunction{iip,false}(f1,f2,g; mass_matrix=I,
                            _func_cache=nothing,analytic=nothing,
                            tgrad = nothing,
                            jac = nothing,
                            jac_prototype = nothing,
                            Wfact = nothing,
                            Wfact_t = nothing,
                            paramjac = nothing,
                            syms = nothing,
                            colorvec = nothing) where iip
  SplitSDEFunction{iip,Any,Any,Any,Any,Any,
                   Any,Any,Any,
                   Any,Any,Any,Any,Any}(
                   f1,f2,g,mass_matrix,_func_cache,analytic,
                   tgrad,jac,jac_prototype,
                   Wfact,Wfact_t,paramjac,syms,colorvec)
end
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

                 if jac === nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                  if iip
                    jac = update_coefficients! #(J,u,p,t)
                  else
                    jac = (u,p,t) -> update_coefficients!(deepcopy(jac_prototype),u,p,t)
                  end
                 end

                 if jac_prototype !== nothing && colorvec === nothing && ArrayInterface.fast_matrix_colors(jac_prototype)
                   _colorvec = ArrayInterface.matrix_colors(jac_prototype)
                 else
                   _colorvec = colorvec
                 end

                 RODEFunction{iip,typeof(f),typeof(mass_matrix),
                 typeof(analytic),typeof(tgrad),
                 typeof(jac),typeof(jac_prototype),typeof(Wfact),typeof(Wfact_t),
                 typeof(paramjac),typeof(syms),typeof(_colorvec)}(
                 f,mass_matrix,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
                 paramjac,syms,_colorvec)
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

                 if jac === nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                  if iip
                    jac = update_coefficients! #(J,u,p,t)
                  else
                    jac = (u,p,t) -> update_coefficients!(deepcopy(jac_prototype),u,p,t)
                  end
                 end

                 if jac_prototype !== nothing && colorvec === nothing && ArrayInterface.fast_matrix_colors(jac_prototype)
                   _colorvec = ArrayInterface.matrix_colors(jac_prototype)
                 else
                   _colorvec = colorvec
                 end

                 RODEFunction{iip,Any,Any,Any,Any,
                 Any,Any,Any,Any,
                 Any,typeof(syms),typeof(_colorvec)}(
                 f,mass_matrix,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
                 paramjac,syms,_colorvec)
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

                 if jac === nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                  if iip
                    jac = update_coefficients! #(J,u,p,t)
                  else
                    jac = (u,p,t) -> update_coefficients!(deepcopy(jac_prototype),u,p,t)
                  end
                 end

                 if jac_prototype !== nothing && colorvec === nothing && ArrayInterface.fast_matrix_colors(jac_prototype)
                   _colorvec = ArrayInterface.matrix_colors(jac_prototype)
                 else
                   _colorvec = colorvec
                 end

                 DAEFunction{iip,typeof(f),typeof(analytic),typeof(tgrad),
                 typeof(jac),typeof(jac_prototype),typeof(Wfact),typeof(Wfact_t),
                 typeof(paramjac),typeof(syms),typeof(_colorvec)}(
                 f,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
                 paramjac,syms,_colorvec)
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

                 if jac === nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                  if iip
                    jac = update_coefficients! #(J,u,p,t)
                  else
                    jac = (u,p,t) -> update_coefficients!(deepcopy(jac_prototype),u,p,t)
                  end
                 end

                 if jac_prototype !== nothing && colorvec === nothing && ArrayInterface.fast_matrix_colors(jac_prototype)
                   _colorvec = ArrayInterface.matrix_colors(jac_prototype)
                 else
                   _colorvec = colorvec
                 end

                 DAEFunction{iip,Any,Any,Any,
                 Any,Any,Any,Any,
                 Any,typeof(syms),typeof(_colorvec)}(
                 f,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
                 paramjac,syms,_colorvec)
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

                 if jac === nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                  if iip
                    jac = update_coefficients! #(J,u,p,t)
                  else
                    jac = (u,p,t) -> update_coefficients!(deepcopy(jac_prototype),u,p,t)
                  end
                 end

                 if jac_prototype !== nothing && colorvec === nothing && ArrayInterface.fast_matrix_colors(jac_prototype)
                   _colorvec = ArrayInterface.matrix_colors(jac_prototype)
                 else
                   _colorvec = colorvec
                 end

                 DDEFunction{iip,typeof(f),typeof(mass_matrix),typeof(analytic),typeof(tgrad),
                 typeof(jac),typeof(jac_prototype),typeof(Wfact),typeof(Wfact_t),
                 typeof(paramjac),typeof(syms),typeof(_colorvec)}(
                 f,mass_matrix,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
                 paramjac,syms,_colorvec)
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

                 if jac === nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
                  if iip
                    jac = update_coefficients! #(J,u,p,t)
                  else
                    jac = (u,p,t) -> update_coefficients!(deepcopy(jac_prototype),u,p,t)
                  end
                 end

                 if jac_prototype !== nothing && colorvec === nothing && ArrayInterface.fast_matrix_colors(jac_prototype)
                   _colorvec = ArrayInterface.matrix_colors(jac_prototype)
                 else
                   _colorvec = colorvec
                 end

                 DDEFunction{iip,Any,Any,Any,Any,
                 Any,Any,Any,Any,
                 Any,typeof(syms),typeof(_colorvec)}(
                 f,mass_matrix,analytic,tgrad,jac,jac_prototype,Wfact,Wfact_t,
                 paramjac,syms,_colorvec)
end
DDEFunction{iip}(f; kwargs...) where iip = DDEFunction{iip,RECOMPILE_BY_DEFAULT}(f; kwargs...)
DDEFunction{iip}(f::DDEFunction; kwargs...) where iip = f
DDEFunction(f; kwargs...) = DDEFunction{isinplace(f, 5),RECOMPILE_BY_DEFAULT}(f; kwargs...)
DDEFunction(f::DDEFunction; kwargs...) = f

########## Existance Functions

# compatibility
has_invW(f::AbstractDiffEqFunction) = false
has_analytic(f::AbstractDiffEqFunction) = __has_analytic(f) && f.analytic !== nothing
has_jac(f::AbstractDiffEqFunction) = __has_jac(f) && f.jac !== nothing
has_tgrad(f::AbstractDiffEqFunction) = __has_tgrad(f) && f.tgrad !== nothing
has_Wfact(f::AbstractDiffEqFunction) = __has_Wfact(f) && f.Wfact !== nothing
has_Wfact_t(f::AbstractDiffEqFunction) = __has_Wfact_t(f) && f.Wfact_t !== nothing
has_paramjac(f::AbstractDiffEqFunction) = __has_paramjac(f) && f.paramjac !== nothing
has_syms(f::AbstractDiffEqFunction) = __has_syms(f) && f.syms !== nothing
has_colorvec(f::AbstractDiffEqFunction) = __has_colorvec(f) && f.colorvec !== nothing

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
