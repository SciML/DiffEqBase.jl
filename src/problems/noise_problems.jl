struct NoiseProblem{N<:AbstractNoiseProcess,T} <: AbstractNoiseProblem
  noise::N
  tspan::T
  seed::UInt64

  @add_kwonly function NoiseProblem(noise,tspan;seed=UInt64(0))
    _tspan = promote_tspan(tspan)
    new{typeof(noise),typeof(_tspan)}(noise,_tspan,seed)
  end
  NoiseProblem{iip}(args...;kwargs...) where {iip} =  NoiseProblem(args...;kwargs...)
end
