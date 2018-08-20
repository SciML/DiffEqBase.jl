struct StandardSDEProblem end

struct SDEProblem{uType,tType,isinplace,P,NP,F,G,C,ND} <: AbstractSDEProblem{uType,tType,isinplace,ND}
  f::F
  g::G
  u0::uType
  tspan::tType
  p::P
  noise::NP
  callback::C
  noise_rate_prototype::ND
  seed::UInt64
  @add_kwonly function SDEProblem{iip}(f::AbstractSDEFunction{iip},g,u0,
          tspan,p=nothing;
          noise_rate_prototype = nothing,
          noise= nothing, seed = UInt64(0),
          callback = nothing) where {iip}
    _tspan = promote_tspan(tspan)

    new{typeof(u0),typeof(_tspan),
        isinplace(f),typeof(p),
        typeof(noise),typeof(f),typeof(f.g),
        typeof(callback),
        typeof(noise_rate_prototype)}(
        f,f.g,u0,_tspan,p,
        noise,callback,
        noise_rate_prototype,seed)
  end

  function SDEProblem{iip}(f,g,u0,tspan,p=nothing;kwargs...) where {iip}
    SDEProblem(convert(SDEFunction{iip},f,g),g,u0,tspan,p;kwargs...)
  end
end

#=
function SDEProblem(f::AbstractSDEFunction,u0,tspan,p=nothing;kwargs...)
  SDEProblem(f,f.g,u0,tspan,p;kwargs...)
end
=#

function SDEProblem(f::AbstractSDEFunction,g,u0,tspan,p=nothing;kwargs...)
  SDEProblem{isinplace(f)}(f,g,u0,tspan,p;kwargs...)
end

function SDEProblem(f,g,u0,tspan,p=nothing;kwargs...)
  SDEProblem(convert(SDEFunction,f,g),g,u0,tspan,p;kwargs...)
end

abstract type AbstractSplitSDEProblem end
struct SplitSDEProblem{iip} <: AbstractSplitSDEProblem end
# u' = Au + f
function SplitSDEProblem(f1,f2,g,u0,tspan,p=nothing;kwargs...)
  SplitSDEProblem(SplitSDEFunction(f1,f2,g),g,u0,tspan,p;kwargs...)
end

SplitSDEProblem(f::SplitSDEFunction,g,u0,tspan,p=nothing;kwargs...) =
  SplitSDEProblem{isinplace(f)}(f,g,u0,tspan,p;kwargs...)

function SplitSDEProblem{iip}(f1,f2,g,u0,tspan,p=nothing;kwargs...) where iip
  SplitSDEProblem(SplitSDEFunction(f1,f2,g),g,u0,tspan,p;kwargs...)
end
function SplitSDEProblem{iip}(f::SplitSDEFunction,g,u0,tspan,p=nothing;
                                     func_cache=nothing,kwargs...) where iip
  if f.cache == nothing && iip
    cache = similar(u0)
    _f = SplitSDEFunction{iip}(f.f1, f.f2, f.g; mass_matrix=f.mass_matrix,
                              _func_cache=cache, analytic=f.analytic)
  else
    _f = f
  end
  SDEProblem(_f,g,u0,tspan,p;kwargs...)
end
