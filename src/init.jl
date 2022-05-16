value(x) = x
promote_u0(u0,p,t0) = u0
promote_tspan(u0,p,tspan,prob,kwargs) = tspan
get_tmp(x) = nothing
isdistribution(u0) = false

function SciMLBase.tmap(args...)
  error("Zygote must be added to differentiate Zygote? If you see this error, report it.")
end

const IS_OPENBLAS = Ref(true)

function __init__()
  @require Measurements="eff96d63-e80a-5855-80a2-b1b0885c5ab7" begin

    promote_u0(u0::AbstractArray{<:Measurements.Measurement},p::AbstractArray{<:Measurements.Measurement},t0) = u0
    promote_u0(u0,p::AbstractArray{<:Measurements.Measurement},t0) = eltype(p).(u0)

    value(x::Type{Measurements.Measurement{T}}) where {T} = T
    value(x::Measurements.Measurement) = Measurements.value(x)

    @inline fastpow(x::Measurements.Measurement, y::Measurements.Measurement) = x^y

    # Support adaptive steps should be errorless
    @inline function ODE_DEFAULT_NORM(u::AbstractArray{<:Measurements.Measurement,N},t) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(t))) / length(u))
    end
    @inline function ODE_DEFAULT_NORM(u::Array{<:Measurements.Measurement,N},t) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(t))) / length(u))
    end
    @inline ODE_DEFAULT_NORM(u::Measurements.Measurement,t) = abs(Measurements.value(u))
  end

  @require MonteCarloMeasurements="0987c9cc-fe09-11e8-30f0-b96dd679fdca" begin

    promote_u0(u0::AbstractArray{<:MonteCarloMeasurements.AbstractParticles},p::AbstractArray{<:MonteCarloMeasurements.AbstractParticles},t0) = u0
    promote_u0(u0,p::AbstractArray{<:MonteCarloMeasurements.AbstractParticles},t0) = eltype(p).(u0)

    value(x::Type{MonteCarloMeasurements.AbstractParticles{T,N}}) where {T,N} = T
    value(x::MonteCarloMeasurements.AbstractParticles) = mean(x.particles)

    @inline fastpow(x::MonteCarloMeasurements.AbstractParticles, y::MonteCarloMeasurements.AbstractParticles) = x^y

    # Support adaptive steps should be errorless
    @inline function ODE_DEFAULT_NORM(u::AbstractArray{<:MonteCarloMeasurements.AbstractParticles,N},t) where {N}
      sqrt(mean(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(t))))
    end
    @inline function ODE_DEFAULT_NORM(u::AbstractArray{<:MonteCarloMeasurements.AbstractParticles,N},t::AbstractArray{<:MonteCarloMeasurements.AbstractParticles,N}) where {N}
      sqrt(mean(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(value.(t)))))
    end
    @inline ODE_DEFAULT_NORM(u::MonteCarloMeasurements.AbstractParticles,t) = abs(value(u))
  end

  @require Unitful="1986cc42-f94f-5a68-af5c-568840ba703d" begin
    # Support adaptive errors should be errorless for exponentiation
    value(x::Type{Unitful.AbstractQuantity{T,D,U}}) where {T,D,U} = T
    value(x::Unitful.AbstractQuantity) = x.val
    @inline function ODE_DEFAULT_NORM(u::AbstractArray{<:Unitful.AbstractQuantity,N},t) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(t))) / length(u))
    end
    @inline function ODE_DEFAULT_NORM(u::Array{<:Unitful.AbstractQuantity,N},t) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(t))) / length(u))
    end
    @inline ODE_DEFAULT_NORM(u::Unitful.AbstractQuantity,t) = abs(value(u))
    @inline UNITLESS_ABS2(x::Unitful.AbstractQuantity) = real(abs2(x)/oneunit(x)*oneunit(x))
  end

  @require GeneralizedGenerated="6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb" begin
    SciMLBase.numargs(::GeneralizedGenerated.RuntimeFn{Args}) where Args = GeneralizedGenerated.from_type(Args) |> length
  end
end
