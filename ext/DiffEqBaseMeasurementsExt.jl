module DiffEqBaseMeasurementsExt

using DiffEqBase
import SciMLBase: value
using Measurements

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
