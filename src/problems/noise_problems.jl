struct NoiseProblem{N<:AbstractNoiseProcess,T} <: AbstractNoiseProblem
  noise::N
  tspan::Tuple{T,T}
  seed::UInt64
end

@add_kwonly NoiseProblem(noise,tspan;seed=UInt64(0)) = NoiseProblem{typeof(noise),
             promote_type(map(typeof,tspan)...)}(
             noise,tspan,seed)
