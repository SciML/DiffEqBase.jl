function __solve end
function __init end
function solve! end

function init(prob::DEProblem,args...;kwargs...)
  _prob = get_concrete_problem(prob,kwargs)
  if haskey(kwargs,:alg)
    __init(_prob,kwargs[:alg],args...;kwargs...)
  else
    __init(_prob,args...;kwargs...)
  end
end

function solve(prob::DEProblem,args...;kwargs...)
  _prob = get_concrete_problem(prob,kwargs)
  if haskey(kwargs,:alg)
    __solve(_prob,kwargs[:alg],args...;kwargs...)
  else
    __solve(_prob,args...;kwargs...)
  end
end

function get_concrete_problem(prob::AbstractSteadyStateProblem,kwargs)
  if typeof(prob.u0) <: Function
    _u0 = prob.u0(prob.p,prob.tspan[1])
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
    _u0 = prob.u0(prob.p,prob.tspan[1])
  else
    _u0 = prob.u0
  end

  __u0 = handle_distribution_u0(_u0)

  remake(prob;u0=__u0,tspan=_tspan)
end

handle_distribution_u0(_u0) = _u0
@require Distributions="31c24e10-a181-5473-b8eb-7969acd0382f" begin
  handle_distribution_u0(_u0::Sampleable) = rand(_u0)
end
