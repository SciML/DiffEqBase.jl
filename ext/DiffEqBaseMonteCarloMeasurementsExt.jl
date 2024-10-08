module DiffEqBaseMonteCarloMeasurementsExt

if isdefined(Base, :get_extension)
    using DiffEqBase
    import DiffEqBase: value
    using MonteCarloMeasurements
else
    using ..DiffEqBase
    import ..DiffEqBase: value
    using ..MonteCarloMeasurements
end

function DiffEqBase.promote_u0(
        u0::AbstractArray{
            <:MonteCarloMeasurements.AbstractParticles,
        },
        p::AbstractArray{<:MonteCarloMeasurements.AbstractParticles},
        t0)
    u0
end
function DiffEqBase.promote_u0(u0,
        p::AbstractArray{<:MonteCarloMeasurements.AbstractParticles},
        t0)
    eltype(p).(u0)
end

DiffEqBase.value(x::Type{MonteCarloMeasurements.AbstractParticles{T, N}}) where {T, N} = T
DiffEqBase.value(x::MonteCarloMeasurements.AbstractParticles) = mean(x.particles)

# Support adaptive steps should be errorless
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::AbstractArray{
            <:MonteCarloMeasurements.AbstractParticles,
            N}, t) where {N}
    sqrt(mean(x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]),
        zip((value(x) for x in u), Iterators.repeated(t))))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::AbstractArray{
            <:MonteCarloMeasurements.AbstractParticles,
            N},
        t::AbstractArray{
            <:MonteCarloMeasurements.AbstractParticles,
            N}) where {N}
    sqrt(mean(x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]),
        zip((value(x) for x in u), Iterators.repeated(value.(t)))))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::MonteCarloMeasurements.AbstractParticles, t)
    abs(value(u))
end

end
