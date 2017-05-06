type NoiseProblem{N<:AbstractNoiseProcess,T} <: AbstractNoiseProblem
  noise::N
  tspan::Tuple{T,T}
end
NoiseProblem(noise,tspan) = NoiseProblem{typeof(noise),promote_type(map(typeof,tspan)...)}(noise,tspan)
