value(x) = x
promote_tspan(u0, p, tspan, prob, kwargs) = tspan
isdistribution(u0) = false

function SciMLBase.tmap(args...)
    error("Zygote must be added to differentiate Zygote? If you see this error, report it.")
end

function __init__()
    @require Measurements="eff96d63-e80a-5855-80a2-b1b0885c5ab7" begin
        function promote_u0(u0::AbstractArray{<:Measurements.Measurement},
                            p::AbstractArray{<:Measurements.Measurement}, t0)
            u0
        end
        promote_u0(u0, p::AbstractArray{<:Measurements.Measurement}, t0) = eltype(p).(u0)

        value(x::Type{Measurements.Measurement{T}}) where {T} = T
        value(x::Measurements.Measurement) = Measurements.value(x)

        @inline fastpow(x::Measurements.Measurement, y::Measurements.Measurement) = x^y

        # Support adaptive steps should be errorless
        @inline function ODE_DEFAULT_NORM(u::AbstractArray{<:Measurements.Measurement, N},
                                          t) where {N}
            sqrt(sum(x -> ODE_DEFAULT_NORM(x[1], x[2]),
                     zip((value(x) for x in u), Iterators.repeated(t))) / length(u))
        end
        @inline function ODE_DEFAULT_NORM(u::Array{<:Measurements.Measurement, N},
                                          t) where {N}
            sqrt(sum(x -> ODE_DEFAULT_NORM(x[1], x[2]),
                     zip((value(x) for x in u), Iterators.repeated(t))) / length(u))
        end
        @inline function ODE_DEFAULT_NORM(u::Measurements.Measurement, t)
            abs(Measurements.value(u))
        end
    end

    @require MonteCarloMeasurements="0987c9cc-fe09-11e8-30f0-b96dd679fdca" begin
        function promote_u0(u0::AbstractArray{<:MonteCarloMeasurements.AbstractParticles},
                            p::AbstractArray{<:MonteCarloMeasurements.AbstractParticles},
                            t0)
            u0
        end
        function promote_u0(u0,
                            p::AbstractArray{<:MonteCarloMeasurements.AbstractParticles},
                            t0)
            eltype(p).(u0)
        end

        value(x::Type{MonteCarloMeasurements.AbstractParticles{T, N}}) where {T, N} = T
        value(x::MonteCarloMeasurements.AbstractParticles) = mean(x.particles)

        @inline function fastpow(x::MonteCarloMeasurements.AbstractParticles,
                                 y::MonteCarloMeasurements.AbstractParticles)
            x^y
        end

        # Support adaptive steps should be errorless
        @inline function ODE_DEFAULT_NORM(u::AbstractArray{
                                                           <:MonteCarloMeasurements.AbstractParticles,
                                                           N}, t) where {N}
            sqrt(mean(x -> ODE_DEFAULT_NORM(x[1], x[2]),
                      zip((value(x) for x in u), Iterators.repeated(t))))
        end
        @inline function ODE_DEFAULT_NORM(u::AbstractArray{
                                                           <:MonteCarloMeasurements.AbstractParticles,
                                                           N},
                                          t::AbstractArray{
                                                           <:MonteCarloMeasurements.AbstractParticles,
                                                           N}) where {N}
            sqrt(mean(x -> ODE_DEFAULT_NORM(x[1], x[2]),
                      zip((value(x) for x in u), Iterators.repeated(value.(t)))))
        end
        @inline function ODE_DEFAULT_NORM(u::MonteCarloMeasurements.AbstractParticles, t)
            abs(value(u))
        end
    end

    @require Unitful="1986cc42-f94f-5a68-af5c-568840ba703d" begin
        # Support adaptive errors should be errorless for exponentiation
        value(x::Type{Unitful.AbstractQuantity{T, D, U}}) where {T, D, U} = T
        value(x::Unitful.AbstractQuantity) = x.val
        @inline function ODE_DEFAULT_NORM(u::AbstractArray{<:Unitful.AbstractQuantity, N},
                                          t) where {N}
            sqrt(sum(x -> ODE_DEFAULT_NORM(x[1], x[2]),
                     zip((value(x) for x in u), Iterators.repeated(t))) / length(u))
        end
        @inline function ODE_DEFAULT_NORM(u::Array{<:Unitful.AbstractQuantity, N},
                                          t) where {N}
            sqrt(sum(x -> ODE_DEFAULT_NORM(x[1], x[2]),
                     zip((value(x) for x in u), Iterators.repeated(t))) / length(u))
        end
        @inline ODE_DEFAULT_NORM(u::Unitful.AbstractQuantity, t) = abs(value(u))
        @inline function UNITLESS_ABS2(x::Unitful.AbstractQuantity)
            real(abs2(x) / oneunit(x) * oneunit(x))
        end
    end

    @require GeneralizedGenerated="6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb" begin function SciMLBase.numargs(::GeneralizedGenerated.RuntimeFn{
                                                                                                                                           Args
                                                                                                                                           }) where {
                                                                                                                                                     Args
                                                                                                                                                     }
        GeneralizedGenerated.from_type(Args) |> length
    end end
end
