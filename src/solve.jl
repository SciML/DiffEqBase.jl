function __solve end
function __init end
function solve! end

function init(prob::DEProblem,args...;kwargs...)
  _prob = get_concrete_problem(prob,kwargs)
  if haskey(kwargs,:alg) && (isempty(args) || args[1] === nothing)
    alg = kwargs[:alg]
    isadaptive(alg) && adaptive_warn(_prob.u0,_prob.tspan)
    __init(_prob,alg,args...;kwargs...)
  elseif !isempty(args) && typeof(args[1]) <: DEAlgorithm
    alg = args[1]
    isadaptive(alg) && adaptive_warn(_prob.u0,_prob.tspan)
    __init(_prob,args...;kwargs...)
  else
    __init(_prob,args...;kwargs...)
  end
end

function solve(prob::DEProblem,args...;kwargs...)
  _prob = get_concrete_problem(prob,kwargs)
  if haskey(kwargs,:alg) && (isempty(args) || args[1] === nothing)
    alg = kwargs[:alg]
    isadaptive(alg) && adaptive_warn(_prob.u0,_prob.tspan)
    __solve(_prob,alg,args...;kwargs...)
  elseif !isempty(args) && typeof(args[1]) <: DEAlgorithm
    alg = args[1]
    isadaptive(alg) && adaptive_warn(_prob.u0,_prob.tspan)
    __solve(_prob,args...;kwargs...)
  else
    __solve(_prob,args...;kwargs...)
  end
end

function get_concrete_problem(prob::AbstractSteadyStateProblem,kwargs)
  if typeof(prob.u0) <: Function
    _u0 = prob.u0(prob.p,Inf)
  else
    _u0 = prob.u0
  end

  __u0 = handle_distribution_u0(_u0)

  remake(prob;u0=__u0)
end

function get_concrete_problem(prob,kwargs)
  if typeof(prob.tspan) <: Function
    _tspan = prob.tspan(prob.p)
  elseif prob.tspan == (nothing,nothing)
    if haskey(kwargs,:tspan)
      _tspan = kwargs.tspan
    else
      error("No tspan is set in the problem or chosen in the init/solve call")
    end
  else
    _tspan = prob.tspan
  end

  if typeof(prob.u0) <: Function
    _u0 = prob.u0(prob.p,_tspan[1])
  else
    _u0 = prob.u0
  end

  __u0 = handle_distribution_u0(_u0)

  remake(prob;u0=__u0,tspan=_tspan)
end

handle_distribution_u0(_u0) = _u0
@require Distributions="31c24e10-a181-5473-b8eb-7969acd0382f" begin
  handle_distribution_u0(_u0::Distributions.Sampleable) = rand(_u0)
end

function adaptive_warn(u0,tspan)
  adaptive_integer_warn(tspan)
  dual_number_warn(u0,tspan)
  measurements_warn(u0,tspan)
end

function adaptive_integer_warn(tspan)
  if eltype(tspan) <: Integer
    @warn("Integer time values are incompatible with adaptive integrators. Utilize floating point numbers instead of integers in this case, i.e. (0.0,1.0) instead of (0,1).")
  end
end

dual_number_warn(u0,tspan) = nothing
@require ForwardDiff="f6369f11-7733-5829-9624-2563aa707210" begin
  function dual_number_warn(u0::AbstractArray{<:Dual},tspan::Tuple{T,T}) T<:Number
    if !(T<:Dual)
      @warn("Both the initial condition and time values must be Dual numbers in order to be compatible with Dual number inputs. Change the element type of tspan to match the element type of u0.")
    end
  end
  function dual_number_warn(u0::AbstractArray{<:Number},tspan::Tuple{T,T}) T<:Dual
    if !(eltype(u0)<:Dual)
      @warn("Both the initial condition and time values must be Dual numbers in order to be compatible with Dual number inputs. Change the element type of u0 to match the element type of tspan.")
    end
  end
end

measurements_warn(u0,tspan) = nothing
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
