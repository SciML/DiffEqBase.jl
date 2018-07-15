struct StandardSDEProblem end

struct SDEProblem{uType,tType,isinplace,P,J,NP,F,G,C,MM,ND} <: AbstractSDEProblem{uType,tType,isinplace,ND}
  f::F
  g::G
  u0::uType
  tspan::tType
  p::P
  jac_prototype::J
  noise::NP
  callback::C
  mass_matrix::MM
  noise_rate_prototype::ND
  seed::UInt64
  @add_kwonly function SDEProblem(f::AbstractSDEFunction,g,u0,
          tspan,p=nothing,problem_type=StandardSDEProblem();
          jac_prototype = nothing,
          noise_rate_prototype = nothing,
          noise= nothing, seed = UInt64(0),
          callback = nothing,mass_matrix=I)
    _tspan = promote_tspan(tspan)
    if mass_matrix == I && typeof(f) <: Tuple
      _mm = ((I for i in 1:length(f))...,)
    else
      _mm = mass_matrix
    end

    new{typeof(u0),typeof(_tspan),
        isinplace(f),typeof(p),typeof(jac_prototype),
        typeof(noise),typeof(f),typeof(f.g),
        typeof(callback),typeof(_mm),
        typeof(noise_rate_prototype)}(
        f,f.g,u0,_tspan,p,jac_prototype,
        noise,callback,_mm,
        noise_rate_prototype,seed)
  end

  function SDEProblem{iip}(f,g,u0,tspan,p=nothing;kwargs...) where {iip}
    SDEProblem(convert(SDEFunction{iip},f,g),g,u0,tspan,p;kwargs...)
  end
end

function SDEProblem(f,g,u0,tspan,p=nothing;kwargs...)
  SDEProblem(convert(SDEFunction,f,g),g,u0,tspan,p;kwargs...)
end

abstract type AbstractSplitSDEProblem end
struct SplitSDEProblem{iip} <: AbstractSplitSDEProblem end
# u' = Au + f
function SplitSDEProblem(f1,f2,g,u0,tspan,p=nothing;kwargs...)
  iip = isinplace(f2,4)
  SplitSDEProblem{iip}(f1,f2,g,u0,tspan,p;kwargs...)
end
function SplitSDEProblem{iip}(f1,f2,g,u0,tspan,p=nothing;
                                     func_cache=nothing,kwargs...) where iip
  iip ? _func_cache = similar(u0) : _func_cache = nothing
  SDEProblem{iip}(SplitFunction{iip}(f1,f2,_func_cache),g,u0,tspan,p,SplitSDEProblem{iip}();kwargs...)
end
