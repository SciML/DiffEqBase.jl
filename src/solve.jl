function __solve end
function __init end
function solve! end

NO_TSPAN_PROBS = Union{AbstractLinearProblem, AbstractNonlinearProblem,
                       AbstractQuadratureProblem,
                       AbstractSteadyStateProblem,AbstractJumpProblem}

has_kwargs(_prob::DEProblem) = has_kwargs(typeof(_prob))
Base.@pure has_kwargs(::Type{T}) where T = :kwargs ∈ fieldnames(T)

function init_call(_prob,args...;kwargs...)
  if has_kwargs(_prob)
    __init(_prob,args...;_prob.kwargs...,kwargs...)
  else
    __init(_prob,args...;kwargs...)
  end
end

function init(prob::DEProblem,args...;kwargs...)
  _prob = get_concrete_problem(prob,kwargs)
  if haskey(kwargs,:alg) && (isempty(args) || args[1] === nothing)
    alg = kwargs[:alg]
    isadaptive(alg) &&
    !(typeof(prob) <: NO_TSPAN_PROBS) &&
    adaptive_warn(_prob.u0,_prob.tspan)
    init_call(_prob,alg,args...;kwargs...)
  elseif !isempty(args) && typeof(args[1]) <: DEAlgorithm
    alg = args[1]
    isadaptive(alg) &&
    !(typeof(prob) <: NO_TSPAN_PROBS) &&
    adaptive_warn(_prob.u0,_prob.tspan)
    init_call(_prob,args...;kwargs...)
  else
    init_call(_prob,args...;kwargs...)
  end
end

function solve_call(_prob,args...;merge_callbacks = true, kwargs...)
  if has_kwargs(_prob)
    if merge_callbacks && haskey(_prob.kwargs,:callback) && haskey(kwargs, :callback)
      kwargs_temp = NamedTuple{Base.diff_names(Base._nt_names(
      values(kwargs)), (:callback,))}(values(kwargs))
      callbacks = NamedTuple{(:callback,)}( [DiffEqBase.CallbackSet(_prob.kwargs[:callback], values(kwargs).callback )] )
      kwargs = merge(kwargs_temp, callbacks)
    end
    kwargs = isempty(_prob.kwargs) ? kwargs : merge(values(_prob.kwargs), kwargs)
  end

  T = Core.Compiler.return_type(__solve,Tuple{typeof(_prob),map(typeof, args)...})

  progress = get(kwargs, :progress, false)
  if progress
    logger = default_logger(Logging.current_logger())
    x = maybe_with_logger(logger) do
      __solve(_prob,args...; kwargs...)
    end
    return x#::T
  else
    __solve(_prob,args...; kwargs...)#::T
  end
end

function solve(prob::DEProblem,args...;sensealg=nothing,
               u0 = nothing, p = nothing,kwargs...)
  u0 = u0 !== nothing ? u0 : prob.u0
  p  = p  !== nothing ? p  : prob.p
  solve_up(prob,sensealg,u0,p,args...;kwargs...)
end

function solve_up(prob::DEProblem,sensealg,u0,p,args...;kwargs...)
  _prob = get_concrete_problem(prob;u0=u0,p=p,kwargs...)
  if haskey(kwargs,:alg) && (isempty(args) || args[1] === nothing)
    alg = kwargs[:alg]
    isadaptive(alg) &&
    !(typeof(prob) <: NO_TSPAN_PROBS) &&
    adaptive_warn(_prob.u0,_prob.tspan)
    solve_call(_prob,alg,args...;kwargs...)
  elseif !isempty(args) && typeof(args[1]) <: DEAlgorithm
    alg = args[1]
    isadaptive(alg) &&
    !(typeof(prob) <: NO_TSPAN_PROBS) &&
    adaptive_warn(_prob.u0,_prob.tspan)
    solve_call(_prob,args...;kwargs...)
  elseif isempty(args) # Default algorithm handling
    !(typeof(prob) <: NO_TSPAN_PROBS) &&
    adaptive_warn(_prob.u0,_prob.tspan)
    solve_call(_prob,args...;kwargs...)
  else
    solve_call(_prob,args...;kwargs...)
  end
end

function solve(prob::EnsembleProblem,args...;kwargs...)
  if isempty(args)
    __solve(prob,nothing,args...;kwargs...)
  else
    __solve(prob,args...;kwargs...)
  end
end

function solve(prob::AbstractNoiseProblem,args...; kwargs...)
  __solve(prob,args...;kwargs...)
end

function get_concrete_problem(prob::AbstractJumpProblem; kwargs...)
  prob
end

function get_concrete_problem(prob::AbstractSteadyStateProblem; kwargs...)
  u0 = get_concrete_u0(prob, Inf, kwargs)
  u0 = promote_u0(u0, prob.p, nothing)
  remake(prob; u0 = u0)
end

function get_concrete_problem(prob::AbstractEnsembleProblem; kwargs...)
  prob
end

function DiffEqBase.solve(prob::PDEProblem,alg::DiffEqBase.DEAlgorithm,args...;
                                          kwargs...)
    solve(prob.prob,alg,args...;kwargs...)
end

function discretize end

function get_concrete_problem(prob; kwargs...)
  p = get_concrete_p(prob, kwargs)
  tspan = get_concrete_tspan(prob, kwargs, p)
  u0 = get_concrete_u0(prob, tspan[1], kwargs)
  u0_promote = promote_u0(u0, p, tspan[1])
  tspan_promote = promote_tspan(u0, p, tspan, prob, kwargs)
  if isconcreteu0(prob, tspan[1], kwargs) && typeof(u0_promote) === typeof(u0) &&
                  prob.tspan == tspan && typeof(tspan) === typeof(tspan_promote)
    return prob
  else
    return remake(prob; u0 = u0_promote, p=p, tspan = tspan_promote)
  end
end

function get_concrete_problem(prob::DDEProblem; kwargs...)
  p = get_concrete_p(prob, kwargs)
  tspan = get_concrete_tspan(prob, kwargs, p)

  u0 = get_concrete_u0(prob, tspan[1], kwargs)

  if prob.constant_lags isa Function
    constant_lags = prob.constant_lags(p)
  else
    constant_lags = prob.constant_lags
  end

  u0 = promote_u0(u0, p, tspan[1])
  tspan = promote_tspan(u0, p, tspan, prob, kwargs)

  remake(prob; u0 = u0, tspan = tspan, p=p, constant_lags = constant_lags)
end

function get_concrete_tspan(prob, kwargs, p)
  if prob.tspan isa Function
    tspan = prob.tspan(p)
  elseif haskey(kwargs, :tspan)
      tspan = kwargs[:tspan]
  elseif prob.tspan === (nothing, nothing)
    error("No tspan is set in the problem or chosen in the init/solve call")
  else
    tspan = prob.tspan
  end

  tspan
end

function isconcreteu0(prob, t0, kwargs)
  !eval_u0(prob.u0) && prob.u0 !== nothing && !isdistribution(prob.u0)
end

function get_concrete_u0(prob, t0, kwargs)
  if eval_u0(prob.u0)
    u0 = prob.u0(prob.p, t0)
  elseif haskey(kwargs,:u0)
    u0 = kwargs[:u0]
  else
    u0 = prob.u0
  end

  handle_distribution_u0(u0)
end

function get_concrete_p(prob, kwargs)
  if haskey(kwargs,:p)
    p = kwargs[:p]
  else
    p = prob.p
  end
end

handle_distribution_u0(_u0) = _u0
eval_u0(u0::Function) = true
eval_u0(u0) = false

"""
$(SIGNATURES)

Check whether the values of `u0` and `tspan` are appropriate for use with
adaptive integrators and emit specific warnings if they are not.
"""
function adaptive_warn(u0,tspan)
  adaptive_integer_warn(tspan)
end

"""
$(SIGNATURES)

Emit a warning about incompatibility with adaptive integers if `tspan` contains
integers.
"""
function adaptive_integer_warn(tspan)
  if eltype(tspan) <: Integer
    @warn("Integer time values are incompatible with adaptive integrators. Utilize floating point numbers instead of integers in this case, i.e. (0.0,1.0) instead of (0,1).")
  end
end

function __solve(prob::DEProblem,args...;default_set=false,second_time=false,kwargs...)
  if second_time
    error("Default algorithm choices require DifferentialEquations.jl. Please specify an algorithm or import DifferentialEquations directly.")
  elseif length(args) > 0 && !(typeof(args[1]) <: Union{Nothing,DEAlgorithm})
    error("Inappropiate solve command. The arguments do not make sense. Likely, you gave an algorithm which does not actually exist (or does not `<:DiffEqBase.DEAlgorithm`)")
  else
    __solve(prob::DEProblem,nothing,args...;default_set=false,second_time=true,kwargs...)
  end
end

################### Differentiation

struct SensitivityADPassThrough <: DiffEqBase.DEAlgorithm end

ZygoteRules.@adjoint function solve_up(prob,sensealg::Union{Nothing,AbstractSensitivityAlgorithm},
                                       u0,p,args...;
                                       kwargs...)
  _solve_adjoint(prob,sensealg,u0,p,args...;kwargs...)
end

function ChainRulesCore.frule(::typeof(solve_up),prob,
                              sensealg::Union{Nothing,AbstractSensitivityAlgorithm},
                              u0,p,args...;
                              kwargs...)
  _solve_forward(prob,sensealg,u0,p,args...;kwargs...)
end

function ChainRulesCore.rrule(::typeof(solve_up),prob,
                              sensealg::Union{Nothing,AbstractSensitivityAlgorithm},
                              u0,p,args...;
                              kwargs...)
  _solve_adjoint(prob,sensealg,u0,p,args...;kwargs...)
end

###
### Legacy Dispatches to be Non-Breaking
###

@deprecate concrete_solve(prob::DiffEqBase.DEProblem,alg::Union{DiffEqBase.DEAlgorithm,Nothing},
                        u0=prob.u0,p=prob.p,args...;kwargs...) solve(prob,alg,args...;u0=u0,p=p,kwargs...)

ZygoteRules.@adjoint function concrete_solve(prob::DiffEqBase.DEProblem,
                              alg::Union{DiffEqBase.DEAlgorithm,Nothing},
                              u0=prob.u0,p=prob.p,args...;
                              sensealg=nothing,
                              kwargs...)
  _concrete_solve_adjoint(prob,alg,sensealg,u0,p,args...;kwargs...)
end

function _solve_adjoint(prob,sensealg,u0,p,args...;kwargs...)
  if isempty(args)
    _concrete_solve_adjoint(prob,nothing,sensealg,u0,p;kwargs...)
  else
    _concrete_solve_adjoint(prob,args[1],sensealg,u0,p,Base.tail(args)...;kwargs...)
  end
end

function _solve_forward(prob,sensealg,u0,p,args...;kwargs...)
  if isempty(args)
    _concrete_solve_forward(prob,nothing,sensealg,u0,p;kwargs...)
  else
    _concrete_solve_forward(prob,args[1],sensealg,u0,p,Base.tail(args)...;kwargs...)
  end
end

function _concrete_solve_adjoint(args...;kwargs...)
  error("No adjoint rules exist. Check that you added `using DiffEqSensitivity`")
end

function _concrete_solve_forward(args...;kwargs...)
  error("No sensitivity rules exist. Check that you added `using DiffEqSensitivity`")
end
