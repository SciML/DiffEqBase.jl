"""
$(TYPEDEF)

TODO
"""
struct NoiseProblem{N<:AbstractNoiseProcess,T,K} <: AbstractNoiseProblem
  noise::N
  tspan::T
  seed::UInt64
  kwargs::K
end

@add_kwonly function NoiseProblem(noise,tspan;seed=UInt64(0),kwargs...)
  _tspan = promote_tspan(tspan)
  NoiseProblem{typeof(noise),typeof(_tspan),typeof(kwargs)}(noise,_tspan,seed,kwargs)
end
