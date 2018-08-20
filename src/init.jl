function __init__()
  @require Distributions="31c24e10-a181-5473-b8eb-7969acd0382f" begin
    handle_distribution_u0(_u0::Distributions.Sampleable) = rand(_u0)
  end

  @require ForwardDiff="f6369f11-7733-5829-9624-2563aa707210" begin
    function dual_number_warn(u0::AbstractArray{<:ForwardDiff.Dual},tspan::Tuple{T,T}) T<:Number
      if !(T<:ForwardDiff.Dual)
        @warn("Both the initial condition and time values must be Dual numbers in order to be compatible with Dual number inputs. Change the element type of tspan to match the element type of u0.")
      end
    end
    function dual_number_warn(u0::AbstractArray{<:Number},tspan::Tuple{T,T}) T<:ForwardDiff.Dual
      if !(eltype(u0)<:ForwardDiff.Dual)
        @warn("Both the initial condition and time values must be Dual numbers in order to be compatible with Dual number inputs. Change the element type of u0 to match the element type of tspan.")
      end
    end
  end

  @require Measurements="eff96d63-e80a-5855-80a2-b1b0885c5ab7" begin
    function measurements_warn(u0::AbstractArray{<:Measurement},tspan::Tuple{T,T}) T<:Number
      if !(T<:Measurement)
        @warn("Both the initial condition and time values must be Dual numbers in order to be compatible with Dual number inputs. Change the element type of tspan to match the element type of u0.")
      end
    end
    function measurements_warn(u0::AbstractArray{<:Number},tspan::Tuple{T,T}) T<:Measurement
      if !(eltype(u0)<:Measurement)
        @warn("Both the initial condition and time values must be Dual numbers in order to be compatible with Dual number inputs. Change the element type of u0 to match the element type of tspan.")
      end
    end
  end
end
