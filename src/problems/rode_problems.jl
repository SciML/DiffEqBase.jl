mutable struct RODEProblem{uType,tType,isinplace,P,NP,F,C,MM,ND} <: AbstractRODEProblem{uType,tType,isinplace,ND}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  p::P
  noise::NP
  callback::C
  mass_matrix::MM
  rand_prototype::ND
  seed::UInt64
  function RODEProblem{iip}(f,u0,tspan,p=nothing;
                       rand_prototype = nothing,
                       noise= nothing, seed = UInt64(0),
                       callback=nothing,mass_matrix=I) where {iip}
  new{typeof(u0),promote_type(map(typeof,tspan)...),
              iip,typeof(p),typeof(noise),typeof(f),typeof(callback),
              typeof(mass_matrix),typeof(rand_prototype)}(
              f,u0,tspan,p,noise,callback,mass_matrix,rand_prototype,seed)
  end
end

function RODEProblem(f,u0,tspan,p=nothing;kwargs...)
  iip = typeof(f)<: Tuple ? isinplace(f[2],5) : isinplace(f,5)
  RODEProblem{iip}(f,u0,tspan,p;kwargs...)
end
