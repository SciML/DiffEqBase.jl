value(x) = x

function __init__()
  @require Distributions="31c24e10-a181-5473-b8eb-7969acd0382f" begin
    handle_distribution_u0(_u0::Distributions.Sampleable) = rand(_u0)
  end

  @require ForwardDiff="f6369f11-7733-5829-9624-2563aa707210" begin

    value(x::Type{ForwardDiff.Dual{T,V,N}}) where {T,V,N} = V
    value(x::ForwardDiff.Dual) = value(ForwardDiff.value(x))

    # Support adaptive with non-dual time
    @inline function ODE_DEFAULT_NORM(u::AbstractArray{<:ForwardDiff.Dual,N},t) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(t))) / length(u))
    end
    @inline function ODE_DEFAULT_NORM(u::Array{<:ForwardDiff.Dual,N},t) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(t))) / length(u))
    end
    @inline ODE_DEFAULT_NORM(u::ForwardDiff.Dual,t) = abs(value(u))

    # When time is dual, it shouldn't drop the duals for adaptivity
    @inline function ODE_DEFAULT_NORM(u::AbstractArray{<:ForwardDiff.Dual,N},t::ForwardDiff.Dual) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((x for x in u),Iterators.repeated(t))) / length(u))
    end
    @inline function ODE_DEFAULT_NORM(u::Array{<:ForwardDiff.Dual,N},t::ForwardDiff.Dual) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((x for x in u),Iterators.repeated(t))) / length(u))
    end
    @inline ODE_DEFAULT_NORM(u::ForwardDiff.Dual,t::ForwardDiff.Dual) = abs(u)

    # Type piracy. Should upstream
    Base.nextfloat(d::ForwardDiff.Dual{T,V,N}) where {T,V,N} = ForwardDiff.Dual{T}(nextfloat(d.value), d.partials)
    Base.prevfloat(d::ForwardDiff.Dual{T,V,N}) where {T,V,N} = ForwardDiff.Dual{T}(prevfloat(d.value), d.partials)
  end

  @require Measurements="eff96d63-e80a-5855-80a2-b1b0885c5ab7" begin

    value(x::Type{Measurements.Measurement{T}}) where {T} = T
    value(x::Measurements.Measurement) = Measurements.value(x)

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

    value(x::Type{MonteCarloMeasurements.AbstractParticles{T,N}}) where {T,N} = T
    value(x::MonteCarloMeasurements.AbstractParticles) = mean(x)

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
  end

  @require Flux="587475ba-b771-5e3f-ad9e-33799f191a9c" begin

    value(x::Type{Flux.Tracker.TrackedReal{T}}) where T = T
    value(x::Type{Flux.Tracker.TrackedArray{T,N,A}}) where {T,N,A} = Array{T,N}
    value(x::Flux.Tracker.TrackedReal)  = x.data
    value(x::Flux.Tracker.TrackedArray) = x.data

    # Support adaptive with non-tracked time
    @inline function ODE_DEFAULT_NORM(u::Flux.Tracker.TrackedArray,t) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(t))) / length(u))
    end
    @inline function ODE_DEFAULT_NORM(u::AbstractArray{<:Flux.Tracker.TrackedReal,N},t) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(t))) / length(u))
    end
    @inline function ODE_DEFAULT_NORM(u::Array{<:Flux.Tracker.TrackedReal,N},t) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(t))) / length(u))
    end
    @inline ODE_DEFAULT_NORM(u::Flux.Tracker.TrackedReal,t) = abs(value(u))

    # Support TrackedReal time, don't drop tracking on the adaptivity there
    @inline function ODE_DEFAULT_NORM(u::Flux.Tracker.TrackedArray,t::Flux.Tracker.TrackedReal) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip(u,Iterators.repeated(t))) / length(u))
    end
    @inline function ODE_DEFAULT_NORM(u::AbstractArray{<:Flux.Tracker.TrackedReal,N},t::Flux.Tracker.TrackedReal) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip(u,Iterators.repeated(t))) / length(u))
    end
    @inline function ODE_DEFAULT_NORM(u::Array{<:Flux.Tracker.TrackedReal,N},t::Flux.Tracker.TrackedReal) where {N}
      sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip(u,Iterators.repeated(t))) / length(u))
    end
    @inline ODE_DEFAULT_NORM(u::Flux.Tracker.TrackedReal,t::Flux.Tracker.TrackedReal) = abs(u)
  end

  # Piracy, should get upstreamed
  @require CuArrays="3a865a2d-5b23-5a0f-bc46-62713ec82fae" begin
    function ldiv!(x::CuArrays.CuArray,_qr::CuArrays.CUSOLVER.CuQR,b::CuArrays.CuArray)
      _x = UpperTriangular(_qr.R) \ (_qr.Q' * reshape(b,length(b),1))
      x .= vec(_x)
    end
  end

end
