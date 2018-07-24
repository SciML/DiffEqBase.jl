mutable struct RODEProblem{uType,tType,isinplace,P,NP,F,C,ND} <: AbstractRODEProblem{uType,tType,isinplace,ND}
  f::F
  u0::uType
  tspan::tType
  p::P
  noise::NP
  callback::C
  rand_prototype::ND
  seed::UInt64
  @add_kwonly function RODEProblem(f::RODEFunction,u0,tspan,p=nothing;
                       rand_prototype = nothing,
                       noise= nothing, seed = UInt64(0),
                       callback=nothing)
  _tspan = promote_tspan(tspan)
  new{typeof(u0),typeof(_tspan),
              isinplace(f),typeof(p),
              typeof(noise),typeof(f),typeof(callback),
              typeof(rand_prototype)}(
              f,u0,_tspan,p,noise,callback,
              rand_prototype,seed)
  end
  function RODEProblem{iip}(f,u0,tspan,p=nothing;kwargs...) where {iip}
    RODEProblem(convert(RODEFunction{iip},f),u0,tspan,p;kwargs...)
  end
end

function RODEProblem(f,u0,tspan,p=nothing;kwargs...)
  RODEProblem(convert(RODEFunction,f),u0,tspan,p;kwargs...)
end
