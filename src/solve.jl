struct EvalFunc{F} <: Function
  f::F
end
(f::EvalFunc)(args...) = f.f(args...)

NO_TSPAN_PROBS = Union{AbstractLinearProblem,AbstractNonlinearProblem,
  AbstractQuadratureProblem,
  AbstractSteadyStateProblem,AbstractJumpProblem}

has_kwargs(_prob::DEProblem) = has_kwargs(typeof(_prob))
Base.@pure has_kwargs(::Type{T}) where {T} = :kwargs ∈ fieldnames(T)

const allowedkeywords = (
  :dense,
  :saveat,
  :save_idxs,
  :tstops,
  :d_discontinuities,
  :save_everystep,
  :save_on,
  :save_start,
  :save_end,
  :initialize_save,
  :adaptive,
  :abstol,
  :reltol,
  :dt,
  :dtmax,
  :dtmin,
  :force_dtmin,
  :internalnorm,
  :controller,
  :gamma,
  :beta1,
  :beta2,
  :qmax,
  :qmin,
  :qsteady_min,
  :qsteady_max,
  :qoldinit,
  :failfactor,
  :calck,
  :alias_u0,
  :maxiters,
  :callback,
  :isoutofdomain,
  :unstable_check,
  :verbose,
  :merge_callbacks,
  :progress,
  :progress_steps,
  :progress_name,
  :progress_message,
  :timeseries_errors,
  :dense_errors,
  :calculate_errors,
  :initializealg,
  :alg,
  :save_noise,
  :delta,
  :seed,
  :alg_hints,
  :kwargshandle,
  :trajectories,
  :batch_size,
  :sensealg,
  :advance_to_tstop,
  :stop_at_next_tstop,
  # These two are from the default algorithm handling
  :default_set,
  :second_time
)

const KWARGWARN_MESSAGE = 
"""
Unrecognized keyword arguments found. Future versions will error.
The only allowed keyword arguments to `solve` are: 
$allowedkeywords
See https://diffeq.sciml.ai/stable/basics/common_solver_opts/ for more details.

Set kwargshandle=KeywordArgError for an error message and more details.
Set kwargshandle=KeywordArgSilent to ignore this message.
"""

const KWARGERROR_MESSAGE = 
"""
Unrecognized keyword arguments found.
The only allowed keyword arguments to `solve` are: 
$allowedkeywords
See https://diffeq.sciml.ai/stable/basics/common_solver_opts/ for more details.
"""

struct CommonKwargError <: Exception 
  kwargs
end

function Base.showerror(io::IO, e::CommonKwargError)
  println(io, KWARGERROR_MESSAGE)
  notin = collect(map(x -> x ∉ allowedkeywords, keys(e.kwargs)))
  unrecognized = collect(keys(e.kwargs))[notin]
  print(io, "Unrecognized keyword arguments: $unrecognized")
end

@enum KeywordArgError KeywordArgWarn KeywordArgSilent

const INCOMPATIBLE_U0_MESSAGE = 
"""
Initial condition incompatible with functional form.
Detected an in-place function with an initial condition of type Number or SArray.
This is incompatible because Numbers cannot be mutated, i.e.
`x = 2.0; y = 2.0; x .= y` will error.

If using a immutable initial condition type, please use the out-of-place form.
I.e. define the function `du=f(u,p,t)` instead of attempting to "mutate" the immutable `du`.

If your differential equation function was defined with multiple dispatches and one is
in-place, then the automatic detection will choose in-place. In this case, override the
choice in the problem constructor, i.e. `ODEProblem{false}(f,u0,tspan,p,kwargs...)`.

For a longer discussion on mutability vs immutability and in-place vs out-of-place, see:
https://diffeq.sciml.ai/stable/tutorials/faster_ode_example/#Example-Accelerating-a-Non-Stiff-Equation:-The-Lorenz-Equation
"""

struct IncompatibleInitialConditionError <: Exception end

function Base.showerror(io::IO, e::IncompatibleInitialConditionError)
  print(io, INCOMPATIBLE_U0_MESSAGE)
end

const NO_DEFAULT_ALGORITHM_MESSAGE = 
"""
Default algorithm choices require DifferentialEquations.jl.
Please specify an algorithm (e.g., `solve(prob, Tsit5())` for an ODE)
or import DifferentialEquations directly.

You can find the list of available solvers at https://diffeq.sciml.ai/stable/solvers/ode_solve/
and its associated pages.
"""

struct NoDefaultAlgorithmError <: Exception end

function Base.showerror(io::IO, e::NoDefaultAlgorithmError)
  print(io, NO_DEFAULT_ALGORITHM_MESSAGE)
end

const NO_TSPAN_MESSAGE = 
"""
No tspan is set in the problem or chosen in the init/solve call
"""

struct NoTspanError <: Exception end

function Base.showerror(io::IO, e::NoTspanError)
  print(io, NO_TSPAN_MESSAGE)
end

const NON_SOLVER_MESSAGE = 
"""
The arguments to solve are incorrect.
The second argument must be a solver choice, `solve(prob,alg)`
where `alg` is a `<: DEAlgorithm`, i.e. `Tsti5()`.

Please double check the arguments being sent to the solver.
"""

struct NonSolverError <: Exception end

function Base.showerror(io::IO, e::NonSolverError)
  print(io, NON_SOLVER_MESSAGE)
end

const PROBSOLVER_PAIRING_MESSAGE = 
"""
Incompatible problem+solver pairing.
For example, this can occur if an ODE solver is passed with an SDEProblem.
Solvers are only capable of handling specific problem types. Please double
check that the chosen pairing is capable for handling the given problems.
"""

struct ProblemSolverPairingError <: Exception 
  prob
  alg
end

function Base.showerror(io::IO, e::ProblemSolverPairingError)
  println(io, PROBSOLVER_PAIRING_MESSAGE)
  println(io, "Problem type: $(SciMLBase.__parameterless_type(typeof(e.prob)))")
  println(io, "Solver type: $(SciMLBase.__parameterless_type(typeof(e.alg)))")
  println(io, "Problem types compatible with the chosen solver: $(compatible_problem_types(e.prob,e.alg))")
end

function compatible_problem_types(prob,alg)
  if alg isa AbstractODEAlgorithm
    ODEProblem
  elseif alg isa AbstractSDEAlgorithm
    (SDEProblem, SDDEProblem)
  elseif alg isa AbstractDDEAlgorithm # StochasticDelayDiffEq.jl just uses the SDE alg
    DDEProblem
  elseif alg isa AbstractDAEAlgorithm
    DAEProblem
  elseif alg isa AbstractSteadyStateAlgorithm
    SteadyStateProblem
  end
end

const DIRECT_AUTODIFF_INCOMPATABILITY_MESSAGE = 
"""
Incompatible solver + automatic differentiation pairing.
The chosen automatic differentiation algorithm requires the ability
for compiler transforms on the code which is only possible on pure-Julia
solvers such as those from OrdinaryDiffEq.jl. Direct differentiation methods 
which require this ability include:

- Direct use of ForwardDiff.jl on the solver
- `ForwardDiffSensitivity`, `ReverseDiffAdjoint`, `TrackerAdjoint`, and `ZygoteAdjoint`
  sensealg choices for adjoint differentiation.

Either switch the choice of solver to a pure Julia method, or change the automatic 
differentiation method to one that does not require such transformations. 

For more details on automatic differentiation, adjoint, and sensitivity analysis
of differential equations, see the documentation page:

https://diffeq.sciml.ai/stable/analysis/sensitivity/
"""

struct DirectAutodiffError <: Exception end

function Base.showerror(io::IO, e::DirectAutodiffError)
  println(io, DIRECT_AUTODIFF_INCOMPATABILITY_MESSAGE)
end

const NONCONCRETE_ELTYPE_MESSAGE = 
"""
Non-concrete element type inside of an `Array` detected.
Arrays with non-concrete element types, such as 
`Array{Union{Float32,Float64}}`, are not supported by the 
differential equation solvers. Anyways, this is bad for
performance so you don't want to be doing this!

If this was a mistake, promote the element types to be
all the same. If this was intentional, for example,
using Unitful.jl with different unit values, then use
an array type which has fast broadcast support for
heterogeneous values such as the ArrayPartition
from RecursiveArrayTools.jl. For example:

```julia
using RecursiveArrayTools
x = ArrayPartition([1.0,2.0],[1f0,2f0])
y = ArrayPartition([3.0,4.0],[3f0,4f0])
x .+ y # fast, stable, and usable as u0 into DiffEq!
```

Element type:
"""

struct NonConcreteEltypeError <: Exception 
  eltype
end

function Base.showerror(io::IO, e::NonConcreteEltypeError)
  print(io, NONCONCRETE_ELTYPE_MESSAGE)
  print(io,e.eltype)
end

const GENERIC_NUMBER_TYPE_ERROR_MESSAGE = 
"""
Non-standard number type (i.e. not Float32, Flaot64,
ComplexF32, or ComplexF64) detected as the element type
for the initial condition or time span. These generic 
number types are only compatible with the pure Julia 
solvers which support generic programming, such as 
OrdinaryDiffEq.jl. The chosen solver does not support 
this functionality. Please double check that the initial 
condition and time span types are correct, and check that
the chosen solver was correct.
"""

struct GenericNumberTypeError <: Exception
  alg
  uType
  tType
end

function Base.showerror(io::IO, e::GenericNumberTypeError)
  println(io, GENERIC_NUMBER_TYPE_ERROR_MESSAGE)
  println(io, "Solver: $(e.alg)")
  println(io, "u0 type: $(e.uType)")
  print(io,"Timespan type: $(e.tType)")
end

function init_call(_prob, args...; merge_callbacks=true, kwargshandle=KeywordArgWarn, kwargs...)

  if has_kwargs(_prob)
    if merge_callbacks && haskey(_prob.kwargs, :callback) && haskey(kwargs, :callback)
      kwargs_temp = NamedTuple{Base.diff_names(Base._nt_names(
          values(kwargs)), (:callback,))}(values(kwargs))
      callbacks = NamedTuple{(:callback,)}((DiffEqBase.CallbackSet(_prob.kwargs[:callback], values(kwargs).callback),))
      kwargs = merge(kwargs_temp, callbacks)
    end
    kwargs = isempty(_prob.kwargs) ? kwargs : merge(values(_prob.kwargs), kwargs)
  end

  checkkwargs(kwargshandle; kwargs...)

  if hasfield(typeof(_prob), :f) && hasfield(typeof(_prob.f), :f) && typeof(_prob.f.f) <: EvalFunc
    Base.invokelatest(__init, _prob, args...; kwargs...)#::T
  else
    __init(_prob, args...; kwargs...)#::T
  end
end

function init(prob::DEProblem, args...; kwargs...)
  if haskey(kwargs, :alg) && (isempty(args) || args[1] === nothing)
    alg = kwargs[:alg]
    _prob = get_concrete_problem(prob, isadaptive(alg); kwargs...)
    check_prob_alg_pairing(_prob, alg)
    init_call(_prob, alg, args...; kwargs...)
  elseif !isempty(args) && typeof(args[1]) <: DEAlgorithm
    alg = args[1]
    _prob = get_concrete_problem(prob, isadaptive(alg); kwargs...)
    check_prob_alg_pairing(_prob, alg)
    init_call(_prob, args...; kwargs...)
  else
    _prob = get_concrete_problem(prob, !(typeof(prob) <: DiscreteProblem); kwargs...)
    init_call(_prob, args...; kwargs...)
  end
end

function solve_call(_prob, args...; merge_callbacks=true, kwargshandle=KeywordArgWarn, kwargs...)
  if has_kwargs(_prob)
    if merge_callbacks && haskey(_prob.kwargs, :callback) && haskey(kwargs, :callback)
      kwargs_temp = NamedTuple{Base.diff_names(Base._nt_names(
          values(kwargs)), (:callback,))}(values(kwargs))
      callbacks = NamedTuple{(:callback,)}((DiffEqBase.CallbackSet(_prob.kwargs[:callback], values(kwargs).callback),))
      kwargs = merge(kwargs_temp, callbacks)
    end
    kwargs = isempty(_prob.kwargs) ? kwargs : merge(values(_prob.kwargs), kwargs)
  end

  checkkwargs(kwargshandle; kwargs...)
  isdefined(_prob, :u0) && !isconcretetype(eltype(_prob.u0)) && 
                          throw(NonConcreteEltypeError(eltype(_prob.u0)))

  if hasfield(typeof(_prob), :f) && hasfield(typeof(_prob.f), :f) && typeof(_prob.f.f) <: EvalFunc
    Base.invokelatest(__solve, _prob, args...; kwargs...)#::T
  else
    __solve(_prob, args...; kwargs...)#::T
  end
end

# save_idxs and saveat are here due to https://github.com/FluxML/Zygote.jl/issues/664
function solve(prob::DEProblem, args...; sensealg=nothing,
  u0=nothing, p=nothing, kwargs...)
  u0 = u0 !== nothing ? u0 : prob.u0
  p = p !== nothing ? p : prob.p
  if sensealg === nothing && haskey(prob.kwargs, :sensealg)
    sensealg = prob.kwargs[:sensealg]
  end
  solve_up(prob, sensealg, u0, p, args...; kwargs...)
end

function solve_up(prob::DEProblem, sensealg, u0, p, args...; kwargs...)

  if haskey(kwargs, :alg) && (isempty(args) || args[1] === nothing)
    alg = kwargs[:alg]
    _alg = prepare_alg(alg, u0, p, prob) 
    _prob = get_concrete_problem(prob, isadaptive(_alg); u0=u0, p=p, kwargs...)
    check_prob_alg_pairing(_prob, alg) # use alg for improved inference
    solve_call(_prob, _alg, args...; kwargs...)
  elseif !isempty(args) && typeof(args[1]) <: DEAlgorithm
    alg = args[1]
    _alg = prepare_alg(alg, u0, p, prob)
    _prob = get_concrete_problem(prob, isadaptive(_alg); u0=u0, p=p, kwargs...)
    check_prob_alg_pairing(_prob, alg)
    solve_call(_prob, _alg, Base.tail(args)...; kwargs...)
  elseif isempty(args) # Default algorithm handling
    _prob = get_concrete_problem(prob, !(typeof(prob) <: DiscreteProblem); u0=u0, p=p, kwargs...)
    solve_call(_prob, args...; kwargs...)
  else
    _prob = get_concrete_problem(prob, !(typeof(prob) <: DiscreteProblem); u0=u0, p=p, kwargs...)
    solve_call(_prob, args...; kwargs...)
  end
end

function solve(prob::EnsembleProblem, args...; kwargs...)
  if isempty(args) || length(args) == 1 && typeof(args[1]) <: EnsembleAlgorithm
    __solve(prob, nothing, args...; kwargs...)
  else
    __solve(prob, args...; kwargs...)
  end
end

function solve(prob::AbstractNoiseProblem, args...; kwargs...)
  __solve(prob, args...; kwargs...)
end

function solve(prob::AbstractJumpProblem, args...; kwargs...)
  __solve(prob, args...; kwargs...)
end

function checkkwargs(kwargshandle; kwargs...)
  if any(x -> x ∉ allowedkeywords, keys(kwargs))
    if kwargshandle == KeywordArgError
      throw(CommonKwargError(kwargs))
    elseif kwargshandle == KeywordArgWarn
      @warn KWARGWARN_MESSAGE
    else
      @assert kwargshandle == KeywordArgSilent
    end
  end
end

@non_differentiable checkkwargs(kwargshandle)

function get_concrete_problem(prob::AbstractJumpProblem, isadapt; kwargs...)
  prob
end

function get_concrete_problem(prob::SteadyStateProblem, isadapt; kwargs...)
  u0 = get_concrete_u0(prob, isadapt, Inf, kwargs)
  u0 = promote_u0(u0, prob.p, nothing)
  remake(prob; u0=u0)
end

function get_concrete_problem(prob::NonlinearProblem, isadapt; kwargs...)
  u0 = get_concrete_u0(prob, isadapt, nothing, kwargs)
  u0 = promote_u0(u0, prob.p, nothing)
  remake(prob; u0=u0)
end

function get_concrete_problem(prob::AbstractEnsembleProblem, isadapt; kwargs...)
  prob
end

function solve(prob::PDEProblem, alg::DiffEqBase.DEAlgorithm, args...;
  kwargs...)
  solve(prob.prob, alg, args...; kwargs...)
end

function init(prob::PDEProblem, alg::DiffEqBase.DEAlgorithm, args...;
  kwargs...)
  init(prob.prob, alg, args...; kwargs...)
end

function get_concrete_problem(prob, isadapt; kwargs...)
  p = get_concrete_p(prob, kwargs)
  tspan = get_concrete_tspan(prob, isadapt, kwargs, p)
  u0 = get_concrete_u0(prob, isadapt, tspan[1], kwargs)
  u0_promote = promote_u0(u0, p, tspan[1])
  f_promote = promote_f(prob.f, u0_promote)
  tspan_promote = promote_tspan(u0_promote, p, tspan, prob, kwargs)
  if isconcreteu0(prob, tspan[1], kwargs) && typeof(u0_promote) === typeof(prob.u0) &&
     prob.tspan == tspan && typeof(prob.tspan) === typeof(tspan_promote) &&
     p === prob.p && f_promote === prob.f
    return prob
  else
    return remake(prob; f=f_promote, u0=u0_promote, p=p, tspan=tspan_promote)
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
     prob.tspan == tspan && typeof(prob.tspan) === typeof(tspan_promote) &&
     p === prob.p && f_promote === prob.f
    return prob
  else
    return remake(prob; f=f_promote, du0=du0_promote, u0=u0_promote, p=p, tspan=tspan_promote)
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

  remake(prob; u0=u0, tspan=tspan, p=p, constant_lags=constant_lags)
end

function promote_f(f::F, u0) where {F}
  # Ensure our jacobian will be of the same type as u0
  uElType = u0 === nothing ? Float64 : eltype(u0)
  if isdefined(f, :jac_prototype) && f.jac_prototype isa AbstractArray
    f = @set f.jac_prototype = similar(f.jac_prototype, uElType)
  end
  return f
end

promote_f(f::SplitFunction, u0) = typeof(f.cache) === typeof(u0) && isinplace(f) ? f : remake(f, cache=zero(u0))
prepare_alg(alg, u0, p, f) = alg

function get_concrete_tspan(prob, isadapt, kwargs, p)
  if prob.tspan isa Function
    tspan = prob.tspan(p)
  elseif haskey(kwargs, :tspan)
    tspan = kwargs[:tspan]
  elseif prob.tspan === (nothing, nothing)
    throw(NoTspanError())
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
  elseif haskey(kwargs, :u0)
    u0 = kwargs[:u0]
  else
    u0 = prob.u0
  end

  isadapt && eltype(u0) <: Integer && (u0 = float.(u0))

  _u0 = handle_distribution_u0(u0)

  if isinplace(prob) && (_u0 isa Number || _u0 isa SArray)
    throw(IncompatibleInitialConditionError())
  end

  _u0
end

function get_concrete_du0(prob, isadapt, t0, kwargs)
  if eval_u0(prob.du0)
    du0 = prob.du0(prob.p, t0)
  elseif haskey(kwargs, :du0)
    du0 = kwargs[:du0]
  else
    du0 = prob.du0
  end

  isadapt && eltype(du0) <: Integer && (du0 = float.(du0))

  _du0 = handle_distribution_u0(du0)

  if isinplace(prob) && (_du0 isa Number || _du0 isa SArray)
    throw(IncompatibleInitialConditionError())
  end

  _du0
end

function get_concrete_p(prob, kwargs)
  if haskey(kwargs, :p)
    p = kwargs[:p]
  else
    p = prob.p
  end
end

handle_distribution_u0(_u0) = _u0
handle_distribution_u0(_u0::Distributions.Sampleable) = rand(_u0)
isdistribution(_u0::Distributions.Sampleable) = true

eval_u0(u0::Function) = true
eval_u0(u0) = false

function __solve(prob::DEProblem, args...; default_set=false, second_time=false, kwargs...)
  if second_time
    throw(NoDefaultAlgorithmError())
  elseif length(args) > 0 && !(typeof(args[1]) <: Union{Nothing,DEAlgorithm})
    throw(NonSolverError())
  else
    __solve(prob::DEProblem, nothing, args...; default_set=false, second_time=true, kwargs...)
  end
end

function check_prob_alg_pairing(prob, alg)
  if prob isa ODEProblem && !(alg isa AbstractODEAlgorithm) ||
     prob isa SDEProblem && !(alg isa AbstractSDEAlgorithm) ||
     prob isa SDDEProblem && !(alg isa AbstractSDEAlgorithm) || 
     prob isa DDEProblem && !(alg isa AbstractDDEAlgorithm) ||
     prob isa DAEProblem && !(alg isa AbstractDAEAlgorithm) ||
     prob isa SteadyStateProblem && !(alg isa AbstractSteadyStateAlgorithm)

     throw(ProblemSolverPairingError(prob, alg))
  end

  if isdefined(prob, :u0) && eltype(prob.u0) <: ForwardDiff.Dual && !SciMLBase.isautodifferentiable(alg)
    throw(DirectAutodiffError())
  end

  # Use automatic differentiability support as a proxy for arbitrary number support
  # This can become more granular if a case where arbitrary number support is different
  # From supporting Dual numbers, but there does not seem to be a real case for that.
  # Check for concrete element type so that the non-concrete case throws a better error
  if !SciMLBase.isautodifferentiable(alg) && (isdefined(prob, :u0) 
                                          && Base.isconcretetype(eltype(prob.u0)) &&
                                          !(eltype(prob.u0) <: Union{Float32,Float64,ComplexF32,ComplexF64})) ||
                                          (isdefined(prob, :tspan) 
                                          && Base.isconcretetype(eltype(prob.tspan)) &&
                                          !(eltype(prob.tspan) <: Union{Float32,Float64,ComplexF32,ComplexF64}))
    
    throw(GenericNumberTypeError(alg, isdefined(prob, :u0) ? prob.u0 : nothing,
                                      isdefined(prob, :tspan) ? prob.tspan : nothing))
  end

end

################### Differentiation

struct SensitivityADPassThrough <: SciMLBase.DEAlgorithm end

function ChainRulesCore.frule(::typeof(solve_up), prob,
  sensealg::Union{Nothing,AbstractSensitivityAlgorithm},
  u0, p, args...;
  kwargs...)
  _solve_forward(prob, sensealg, u0, p, args...; kwargs...)
end

function ChainRulesCore.rrule(::typeof(solve_up), prob::SciMLBase.DEProblem,
  sensealg::Union{Nothing,AbstractSensitivityAlgorithm},
  u0, p, args...;
  kwargs...)
  _solve_adjoint(prob, sensealg, u0, p, args...; kwargs...)
end

###
### Legacy Dispatches to be Non-Breaking
###

@deprecate concrete_solve(prob::SciMLBase.DEProblem, alg::Union{SciMLBase.DEAlgorithm,Nothing},
  u0=prob.u0, p=prob.p, args...; kwargs...) solve(prob, alg, args...; u0=u0, p=p, kwargs...)

function _solve_adjoint(prob, sensealg, u0, p, args...; merge_callbacks=true, kwargs...)

  _prob = if haskey(kwargs, :alg) && (isempty(args) || args[1] === nothing)
    alg = kwargs[:alg]
    get_concrete_problem(prob, isadaptive(alg); u0=u0, p=p, kwargs...)
  elseif !isempty(args) && typeof(args[1]) <: DEAlgorithm
    alg = args[1]
    get_concrete_problem(prob, isadaptive(alg); u0=u0, p=p, kwargs...)
  elseif isempty(args) # Default algorithm handling
    get_concrete_problem(prob, !(typeof(prob) <: DiscreteProblem); u0=u0, p=p, kwargs...)
  else
    get_concrete_problem(prob, !(typeof(prob) <: DiscreteProblem); u0=u0, p=p, kwargs...)
  end

  if has_kwargs(_prob)
    if merge_callbacks && haskey(_prob.kwargs, :callback) && haskey(kwargs, :callback)
      kwargs_temp = NamedTuple{Base.diff_names(Base._nt_names(
          values(kwargs)), (:callback,))}(values(kwargs))
      callbacks = NamedTuple{(:callback,)}((DiffEqBase.CallbackSet(_prob.kwargs[:callback], values(kwargs).callback),))
      kwargs = merge(kwargs_temp, callbacks)
    end
    kwargs = isempty(_prob.kwargs) ? kwargs : merge(values(_prob.kwargs), kwargs)
  end

  if isempty(args)
    _concrete_solve_adjoint(_prob, nothing, sensealg, u0, p; kwargs...)
  else
    _concrete_solve_adjoint(_prob, args[1], sensealg, u0, p, Base.tail(args)...; kwargs...)
  end
end

function _solve_forward(prob, sensealg, u0, p, args...; merge_callbacks=true, kwargs...)

  _prob = if haskey(kwargs, :alg) && (isempty(args) || args[1] === nothing)
    alg = kwargs[:alg]
    get_concrete_problem(prob, isadaptive(alg); u0=u0, p=p, kwargs...)
  elseif !isempty(args) && typeof(args[1]) <: DEAlgorithm
    alg = args[1]
    get_concrete_problem(prob, isadaptive(alg); u0=u0, p=p, kwargs...)
  elseif isempty(args) # Default algorithm handling
    get_concrete_problem(prob, !(typeof(prob) <: DiscreteProblem); u0=u0, p=p, kwargs...)
  else
    get_concrete_problem(prob, !(typeof(prob) <: DiscreteProblem); u0=u0, p=p, kwargs...)
  end

  if has_kwargs(_prob)
    if merge_callbacks && haskey(_prob.kwargs, :callback) && haskey(kwargs, :callback)
      kwargs_temp = NamedTuple{Base.diff_names(Base._nt_names(
          values(kwargs)), (:callback,))}(values(kwargs))
      callbacks = NamedTuple{(:callback,)}((DiffEqBase.CallbackSet(_prob.kwargs[:callback], values(kwargs).callback),))
      kwargs = merge(kwargs_temp, callbacks)
    end
    kwargs = isempty(_prob.kwargs) ? kwargs : merge(values(_prob.kwargs), kwargs)
  end

  if isempty(args)
    _concrete_solve_forward(prob, nothing, sensealg, u0, p; kwargs...)
  else
    _concrete_solve_forward(prob, args[1], sensealg, u0, p, Base.tail(args)...; kwargs...)
  end
end

function _concrete_solve_adjoint(args...; kwargs...)
  error("No adjoint rules exist. Check that you added `using DiffEqSensitivity`")
end

function _concrete_solve_forward(args...; kwargs...)
  error("No sensitivity rules exist. Check that you added `using DiffEqSensitivity`")
end
