struct EvalFunc{F} <: Function
  f::F
end
(f::EvalFunc)(args...) = f.f(args...)

NO_TSPAN_PROBS = Union{AbstractLinearProblem, AbstractNonlinearProblem,
                       AbstractQuadratureProblem,
                       AbstractSteadyStateProblem,AbstractJumpProblem}

has_kwargs(_prob::DEProblem) = has_kwargs(typeof(_prob))
Base.@pure has_kwargs(::Type{T}) where T = :kwargs ∈ fieldnames(T)

function init_call(_prob,args...;merge_callbacks = true,kwargs...)

  if has_kwargs(_prob)
    if merge_callbacks && haskey(_prob.kwargs,:callback) && haskey(kwargs, :callback)
      kwargs_temp = NamedTuple{Base.diff_names(Base._nt_names(
      values(kwargs)), (:callback,))}(values(kwargs))
      callbacks = NamedTuple{(:callback,)}( [DiffEqBase.CallbackSet(_prob.kwargs[:callback], values(kwargs).callback )] )
      kwargs = merge(kwargs_temp, callbacks)
    end
    kwargs = isempty(_prob.kwargs) ? kwargs : merge(values(_prob.kwargs), kwargs)
  end

  if hasfield(typeof(_prob),:f) && hasfield(typeof(_prob.f),:f) && typeof(_prob.f.f) <: EvalFunc
    Base.invokelatest(__init,_prob,args...; kwargs...)#::T
  else
    __init(_prob,args...;kwargs...)#::T
  end
end

function init(prob::DEProblem,args...;kwargs...)
  if haskey(kwargs,:alg) && (isempty(args) || args[1] === nothing)
    alg = kwargs[:alg]
    _prob = get_concrete_problem(prob,isadaptive(alg);kwargs...)
    init_call(_prob,alg,args...;kwargs...)
  elseif !isempty(args) && typeof(args[1]) <: DEAlgorithm
    alg = args[1]
    _prob = get_concrete_problem(prob,isadaptive(alg);kwargs...)
    init_call(_prob,args...;kwargs...)
  else
    _prob = get_concrete_problem(prob,!(typeof(prob)<:DiscreteProblem);kwargs...)
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

  if hasfield(typeof(_prob),:f) && hasfield(typeof(_prob.f),:f) && typeof(_prob.f.f) <: EvalFunc
    Base.invokelatest(__solve,_prob,args...; kwargs...)#::T
  else
    __solve(_prob,args...;kwargs...)#::T
  end
end

# save_idxs and saveat are here due to https://github.com/FluxML/Zygote.jl/issues/664
function solve(prob::DEProblem,args...;sensealg=nothing,
               u0 = nothing, p = nothing, kwargs...)
  u0 = u0 !== nothing ? u0 : prob.u0
  p  = p  !== nothing ? p  : prob.p
  solve_up(prob,sensealg,u0,p,args...;kwargs...)
end

function solve_up(prob::DEProblem,sensealg,u0,p,args...;kwargs...)

  if haskey(kwargs,:alg) && (isempty(args) || args[1] === nothing)
    alg = kwargs[:alg]
    _prob = get_concrete_problem(prob,isadaptive(alg);u0=u0,p=p,kwargs...)
    solve_call(_prob,alg,args...;kwargs...)
  elseif !isempty(args) && typeof(args[1]) <: DEAlgorithm
    alg = args[1]
    _prob = get_concrete_problem(prob,isadaptive(alg);u0=u0,p=p,kwargs...)
    solve_call(_prob,args...;kwargs...)
  elseif isempty(args) # Default algorithm handling
    _prob = get_concrete_problem(prob,!(typeof(prob)<:DiscreteProblem);u0=u0,p=p,kwargs...)
    solve_call(_prob,args...;kwargs...)
  else
    _prob = get_concrete_problem(prob,!(typeof(prob)<:DiscreteProblem);u0=u0,p=p,kwargs...)
    solve_call(_prob,args...;kwargs...)
  end
end

function solve(prob::EnsembleProblem,args...;kwargs...)
  if isempty(args) || length(args) == 1 && typeof(args[1])<:EnsembleAlgorithm
    __solve(prob,nothing,args...;kwargs...)
  else
    __solve(prob,args...;kwargs...)
  end
end

function solve(prob::AbstractNoiseProblem,args...; kwargs...)
  __solve(prob,args...;kwargs...)
end

function solve(prob::AbstractJumpProblem,args...; kwargs...)
  __solve(prob,args...;kwargs...)
end

function get_concrete_problem(prob::AbstractJumpProblem, isadapt; kwargs...)
  prob
end

function get_concrete_problem(prob::SteadyStateProblem, isadapt; kwargs...)
  u0 = get_concrete_u0(prob, isadapt, Inf, kwargs)
  u0 = promote_u0(u0, prob.p, nothing)
  remake(prob; u0 = u0)
end

function get_concrete_problem(prob::NonlinearProblem, isadapt; kwargs...)
  u0 = get_concrete_u0(prob, isadapt, nothing, kwargs)
  u0 = promote_u0(u0, prob.p, nothing)
  remake(prob; u0 = u0)
end

function get_concrete_problem(prob::AbstractEnsembleProblem, isadapt; kwargs...)
  prob
end

function solve(prob::PDEProblem,alg::DiffEqBase.DEAlgorithm,args...;
                                          kwargs...)
    solve(prob.prob,alg,args...;kwargs...)
end

function init(prob::PDEProblem,alg::DiffEqBase.DEAlgorithm,args...;
                                          kwargs...)
    init(prob.prob,alg,args...;kwargs...)
end

function get_concrete_problem(prob, isadapt; kwargs...)
  p = get_concrete_p(prob, kwargs)
  tspan = get_concrete_tspan(prob, isadapt, kwargs, p)
  u0 = get_concrete_u0(prob, isadapt, tspan[1], kwargs)
  u0_promote = promote_u0(u0, p, tspan[1])
  f_promote = promote_f(prob.f, u0_promote)
  tspan_promote = promote_tspan(u0_promote, p, tspan, prob, kwargs)
  if isconcreteu0(prob, tspan[1], kwargs) && typeof(u0_promote) === typeof(prob.u0) &&
                  prob.tspan == tspan && typeof(tspan) === typeof(tspan_promote) &&
                  p === prob.p && f_promote === prob.f
    return prob
  else
    return remake(prob; f = f_promote, u0 = u0_promote, p=p, tspan = tspan_promote)
  end
end

function get_concrete_problem(prob::DAEProblem, isadapt; kwargs...)

  p = get_concrete_p(prob, kwargs)
  tspan = get_concrete_tspan(prob, isadapt, kwargs, p)
  u0 = get_concrete_u0(prob, isadapt, tspan[1], kwargs)
  du0 = get_concrete_du0(prob, isadapt, tspan[1], kwargs)

  u0_promote = promote_u0(u0, p, tspan[1])
  du0_promote = promote_u0(du0, p, tspan[1])

  f_promote = promote_f(prob.f, u0_promote)
  tspan_promote = promote_tspan(u0_promote, p, tspan, prob, kwargs)
  if isconcreteu0(prob, tspan[1], kwargs) && typeof(u0_promote) === typeof(prob.u0) &&
     isconcretedu0(prob, tspan[1], kwargs) && typeof(du0_promote) === typeof(prob.du0) &&
                  prob.tspan == tspan && typeof(tspan) === typeof(tspan_promote) &&
                  p === prob.p && f_promote === prob.f
    return prob
  else
    return remake(prob; f = f_promote, du0 = du0_promote, u0 = u0_promote, p=p, tspan = tspan_promote)
  end
end

function get_concrete_problem(prob::DDEProblem, isadapt; kwargs...)
  p = get_concrete_p(prob, kwargs)
  tspan = get_concrete_tspan(prob, isadapt, kwargs, p)
  u0 = get_concrete_u0(prob, isadapt, tspan[1], kwargs)

  if prob.constant_lags isa Function
    constant_lags = prob.constant_lags(p)
  else
    constant_lags = prob.constant_lags
  end

  u0 = promote_u0(u0, p, tspan[1])
  tspan = promote_tspan(u0, p, tspan, prob, kwargs)

  remake(prob; u0 = u0, tspan = tspan, p=p, constant_lags = constant_lags)
end

function promote_f(f,u0)
    # Ensure our jacobian will be of the same type as u0
    uElType = u0 === nothing ? Float64 : eltype(u0)
    if isdefined(f, :jac_prototype) && f.jac_prototype isa AbstractArray
        f = @set f.jac_prototype = similar(f.jac_prototype, uElType)
    end
    return f
end

promote_f(f::SplitFunction,u0) = typeof(f.cache) === typeof(u0) && isinplace(f) ? f : remake(f,cache=zero(u0))

function get_concrete_tspan(prob, isadapt, kwargs, p)
  if prob.tspan isa Function
    tspan = prob.tspan(p)
  elseif haskey(kwargs, :tspan)
      tspan = kwargs[:tspan]
  elseif prob.tspan === (nothing, nothing)
    error("No tspan is set in the problem or chosen in the init/solve call")
  else
    tspan = prob.tspan
  end

  isadapt && eltype(tspan) <: Integer && (tspan = float.(tspan))

  tspan
end

function isconcreteu0(prob, t0, kwargs)
  !eval_u0(prob.u0) && prob.u0 !== nothing && !isdistribution(prob.u0)
end

function isconcretedu0(prob, t0, kwargs)
  !eval_u0(prob.u0) && prob.du0 !== nothing && !isdistribution(prob.du0)
end

function get_concrete_u0(prob, isadapt, t0, kwargs)
  if eval_u0(prob.u0)
    u0 = prob.u0(prob.p, t0)
  elseif haskey(kwargs,:u0)
    u0 = kwargs[:u0]
  else
    u0 = prob.u0
  end

  isadapt && eltype(u0) <: Integer && (u0 = float.(u0))

  handle_distribution_u0(u0)
end

function get_concrete_du0(prob, isadapt, t0, kwargs)
  if eval_u0(prob.du0)
    du0 = prob.du0(prob.p, t0)
  elseif haskey(kwargs,:du0)
    du0 = kwargs[:du0]
  else
    du0 = prob.du0
  end

  isadapt && eltype(du0) <: Integer && (du0 = float.(du0))

  handle_distribution_u0(du0)
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

function __solve(prob::DEProblem,args...;default_set=false,second_time=false,kwargs...)
  if second_time
    error("""
      Default algorithm choices require DifferentialEquations.jl.
      Please specify an algorithm (e.g., `solve(prob, Tsit5())` for an ODE)
      or import DifferentialEquations directly.

      You can find the list of available solvers at https://diffeq.sciml.ai/stable/solvers/ode_solve/
      and its associated pages.
      """
    )
  elseif length(args) > 0 && !(typeof(args[1]) <: Union{Nothing,DEAlgorithm})
    error("Inappropiate solve command. The arguments do not make sense. Likely, you gave an algorithm which does not actually exist (or does not `<:DiffEqBase.DEAlgorithm`)")
  else
    __solve(prob::DEProblem,nothing,args...;default_set=false,second_time=true,kwargs...)
  end
end

################### Differentiation

struct SensitivityADPassThrough <: SciMLBase.DEAlgorithm end

function ChainRulesCore.frule(::typeof(solve_up),prob,
                              sensealg::Union{Nothing,AbstractSensitivityAlgorithm},
                              u0,p,args...;
                              kwargs...)
  _solve_forward(prob,sensealg,u0,p,args...;kwargs...)
end

function ChainRulesCore.rrule(::typeof(solve_up),prob::SciMLBase.DEProblem,
                              sensealg::Union{Nothing,AbstractSensitivityAlgorithm},
                              u0,p,args...;
                              kwargs...)
  _solve_adjoint(prob,sensealg,u0,p,args...;kwargs...)
end

###
### Legacy Dispatches to be Non-Breaking
###

@deprecate concrete_solve(prob::SciMLBase.DEProblem,alg::Union{SciMLBase.DEAlgorithm,Nothing},
                        u0=prob.u0,p=prob.p,args...;kwargs...) solve(prob,alg,args...;u0=u0,p=p,kwargs...)

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
