struct NoiseProblem{N<:AbstractNoiseProcess,T} <: AbstractNoiseProblem
  noise::N
  tspan::T
  seed::UInt64
end

@add_kwonly function NoiseProblem(noise,tspan;seed=UInt64(0))
  _tspan = promote_tspan(tspan)
  NoiseProblem{typeof(noise),typeof(_tspan)}(noise,_tspan,seed)
end
