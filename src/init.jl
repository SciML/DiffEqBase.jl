function __init__()
  @require Distributions="31c24e10-a181-5473-b8eb-7969acd0382f" begin
    handle_distribution_u0(_u0::Distributions.Sampleable) = rand(_u0)
  end

  @require ForwardDiff="f6369f11-7733-5829-9624-2563aa707210" begin
    @inline function ODE_DEFAULT_NORM(u::AbstractArray{<:ForwardDiff.Dual,N}) where {N}
      sqrt(sum(ODE_DEFAULT_NORM,(ForwardDiff.value(x) for x in u)) / length(u))
    end
    @inline function ODE_DEFAULT_NORM(u::Array{<:ForwardDiff.Dual,N}) where {N}
      sqrt(sum(ODE_DEFAULT_NORM,(ForwardDiff.value(x) for x in u)) / length(u))
    end
    @inline ODE_DEFAULT_NORM(u::ForwardDiff.Dual) = abs(ForwardDiff.value(u))
  end

  @require Measurements="eff96d63-e80a-5855-80a2-b1b0885c5ab7" begin
    @inline function ODE_DEFAULT_NORM(u::AbstractArray{<:Measurements.Measurement,N}) where {N}
      sqrt(sum(ODE_DEFAULT_NORM,(Measurements.value(x) for x in u)) / length(u))
    end
    @inline function ODE_DEFAULT_NORM(u::Array{<:Measurements.Measurement,N}) where {N}
      sqrt(sum(ODE_DEFAULT_NORM,(Measurements.value(x) for x in u)) / length(u))
    end
    @inline ODE_DEFAULT_NORM(u::Measurements.Measurement) = abs(Measurements.value(u))
  end

  @require Unitful="1986cc42-f94f-5a68-af5c-568840ba703d" begin
    value(x::Unitful.Quantity) = x.val
    @inline function ODE_DEFAULT_NORM(u::AbstractArray{<:Unitful.Quantity,N}) where {N}
      sqrt(sum(ODE_DEFAULT_NORM,(value(x) for x in u)) / length(u))
    end
    @inline function ODE_DEFAULT_NORM(u::Array{<:Unitful.Quantity,N}) where {N}
      sqrt(sum(ODE_DEFAULT_NORM,(value(x) for x in u)) / length(u))
    end
    @inline ODE_DEFAULT_NORM(u::Unitful.Quantity) = abs(value(u))
  end

  @require Flux="587475ba-b771-5e3f-ad9e-33799f191a9c" begin
    value(x::Flux.Tracker.TrackedReal)  = x.data
    value(x::Flux.Tracker.TrackedArray) = x.data
    @inline function ODE_DEFAULT_NORM(u::Flux.Tracker.TrackedArray) where {N}
      sqrt(sum(ODE_DEFAULT_NORM,(value(x) for x in u)) / length(u))
    end
    @inline ODE_DEFAULT_NORM(u::Flux.Tracker.TrackedReal) = abs(value(u))
  end

end
