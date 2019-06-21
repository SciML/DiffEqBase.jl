function __solve end
function __init end
function solve! end

NO_TSPAN_PROBS = Union{AbstractSteadyStateProblem,AbstractJumpProblem}

function init(prob::DEProblem,args...;kwargs...)
  _prob = get_concrete_problem(prob,kwargs)
  if haskey(kwargs,:alg) && (isempty(args) || args[1] === nothing)
    alg = kwargs[:alg]
    isadaptive(alg) &&
    !(typeof(prob) <: NO_TSPAN_PROBS) &&
    adaptive_warn(_prob.u0,_prob.tspan)
    __init(_prob,alg,args...;kwargs...)
  elseif !isempty(args) && typeof(args[1]) <: DEAlgorithm
    alg = args[1]
    isadaptive(alg) &&
    !(typeof(prob) <: NO_TSPAN_PROBS) &&
    adaptive_warn(_prob.u0,_prob.tspan)
    __init(_prob,args...;kwargs...)
  else
    __init(_prob,args...;kwargs...)
  end
end

function solve(prob::DEProblem,args...;kwargs...)
  _prob = get_concrete_problem(prob,kwargs)
  if haskey(kwargs,:alg) && (isempty(args) || args[1] === nothing)
    alg = kwargs[:alg]
    isadaptive(alg) &&
    !(typeof(prob) <: NO_TSPAN_PROBS) &&
    adaptive_warn(_prob.u0,_prob.tspan)
    __solve(_prob,alg,args...;kwargs...)
  elseif !isempty(args) && typeof(args[1]) <: DEAlgorithm
    alg = args[1]
    isadaptive(alg) &&
    !(typeof(prob) <: NO_TSPAN_PROBS) &&
    adaptive_warn(_prob.u0,_prob.tspan)
    __solve(_prob,args...;kwargs...)
  else
    __solve(_prob,args...;kwargs...)
  end
end

function get_concrete_problem(prob::AbstractJumpProblem,kwargs)
  prob
end

function get_concrete_problem(prob::AbstractSteadyStateProblem, kwargs)
  u0 = get_concrete_u0(prob, Inf)
  remake(prob; u0 = u0)
end

function get_concrete_problem(prob::AbstractMonteCarloProblem, kwargs)
  prob
end

function get_concrete_problem(prob, kwargs)
  tspan = get_concrete_tspan(prob, kwargs)
  u0 = get_concrete_u0(prob, tspan[1])
  remake(prob; u0 = u0, tspan = tspan)
end

function get_concrete_problem(prob::DDEProblem, kwargs)
  tspan = get_concrete_tspan(prob, kwargs)

  u0 = get_concrete_u0(prob, tspan[1])

  if prob.constant_lags isa Function
    constant_lags = prob.constant_lags(prob.p)
  else
    constant_lags = prob.constant_lags
  end

  remake(prob; u0 = u0, tspan = tspan, constant_lags = constant_lags)
end

function get_concrete_tspan(prob, kwargs)
  if prob.tspan isa Function
    tspan = prob.tspan(prob.p)
  elseif prob.tspan == (nothing, nothing)
    if haskey(kwargs, :tspan)
      tspan = kwargs.tspan
    else
      error("No tspan is set in the problem or chosen in the init/solve call")
    end
  else
    tspan = prob.tspan
  end

  tspan
end

function get_concrete_u0(prob, t0)
  if prob.u0 isa Function
    u0 = prob.u0(prob.p, t0)
  else
    u0 = prob.u0
  end

  handle_distribution_u0(u0)
end

handle_distribution_u0(_u0) = _u0

function adaptive_warn(u0,tspan)
  adaptive_integer_warn(tspan)
end

function adaptive_integer_warn(tspan)
  if eltype(tspan) <: Integer
    @warn("Integer time values are incompatible with adaptive integrators. Utilize floating point numbers instead of integers in this case, i.e. (0.0,1.0) instead of (0,1).")
  end
end
