mutable struct RODEProblem{uType,tType,isinplace,P,NP,F,K,ND} <: AbstractRODEProblem{uType,tType,isinplace,ND}
  f::F
  u0::uType
  tspan::tType
  p::P
  noise::NP
  kwargs::K
  rand_prototype::ND
  seed::UInt64
  @add_kwonly function RODEProblem{iip}(f::RODEFunction{iip},u0,tspan,p=NullParameters();
                       rand_prototype = nothing,
                       noise= nothing, seed = UInt64(0),
                       kwargs...) where {iip}
  _tspan = promote_tspan(tspan)
  new{typeof(u0),typeof(_tspan),
              isinplace(f),typeof(p),
              typeof(noise),typeof(f),typeof(kwargs),
              typeof(rand_prototype)}(
              f,u0,_tspan,p,noise,kwargs,
              rand_prototype,seed)
  end
  function RODEProblem{iip}(f,u0,tspan,p=NullParameters();kwargs...) where {iip}
    RODEProblem(convert(RODEFunction{iip},f),u0,tspan,p;kwargs...)
  end
end

function RODEProblem(f::RODEFunction,u0,tspan,p=NullParameters();kwargs...)
  RODEProblem{isinplace(f)}(f,u0,tspan,p;kwargs...)
end

function RODEProblem(f,u0,tspan,p=NullParameters();kwargs...)
  RODEProblem(convert(RODEFunction,f),u0,tspan,p;kwargs...)
end
