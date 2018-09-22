function __init__()
  @require Distributions="31c24e10-a181-5473-b8eb-7969acd0382f" begin
    handle_distribution_u0(_u0::Distributions.Sampleable) = rand(_u0)
  end

  @require ForwardDiff="f6369f11-7733-5829-9624-2563aa707210" begin
    @inline function ODE_DEFAULT_NORM(u::AbstractArray{<:ForwardDiff.Dual,N}) where {N}
      sqrt(sum(ODE_DEFAULT_NORM,(ForwardDiff.value(x) for x in u)) / length(u))
    end
  end

  @require Measurements="eff96d63-e80a-5855-80a2-b1b0885c5ab7" begin
    @inline function ODE_DEFAULT_NORM(u::AbstractArray{<:ForwardDiff.Dual,N}) where {N}
      sqrt(sum(ODE_DEFAULT_NORM,(ForwardDiff.value(x) for x in u)) / length(u))
    end
  end

end
