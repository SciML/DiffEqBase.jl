module DiffEqBaseMeasurementsExt

if isdefined(Base, :get_extension)
    using DiffEqBase
    import DiffEqBase: value
    using Measurements
else
    using ..DiffEqBase
    import ..DiffEqBase: value
    using ..Measurements
end

function DiffEqBase.promote_u0(u0::AbstractArray{<:Measurements.Measurement},
        p::AbstractArray{<:Measurements.Measurement}, t0)
    u0
end
DiffEqBase.promote_u0(u0, p::AbstractArray{<:Measurements.Measurement}, t0) = eltype(p).(u0)

value(x::Type{Measurements.Measurement{T}}) where {T} = T
value(x::Measurements.Measurement) = Measurements.value(x)

# Support adaptive steps should be errorless
@inline function DiffEqBase.ODE_DEFAULT_NORM(
        u::AbstractArray{
            <:Measurements.Measurement,
            N
        },
        t) where {N}
    sqrt(sum(x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]),
        zip((value(x) for x in u), Iterators.repeated(t))) / length(u))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::Array{<:Measurements.Measurement, N},
        t) where {N}
    sqrt(sum(x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]),
        zip((value(x) for x in u), Iterators.repeated(t))) / length(u))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::Measurements.Measurement, t)
    abs(Measurements.value(u))
end

end
