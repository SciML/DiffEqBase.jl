function __solve end
function __init end
function solve! end

function init(prob::DEProblem,args...;kwargs...)
  _prob = get_concrete_problem(prob,kwargs)
  __init(_prob,alg,args...;kwargs...)
end

function solve(prob::DEProblem,args...;kwargs...)
  _prob = get_concrete_problem(prob,kwargs)
  __solve(_prob,alg,args...;kwargs...)
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

  remake(prob;u0=_u0,tspan=_tspan)
end
