struct EvalFunc{F} <: Function
    f::F
end
(f::EvalFunc)(args...) = f.f(args...)

NO_TSPAN_PROBS = Union{AbstractLinearProblem, AbstractNonlinearProblem,
    AbstractIntegralProblem, AbstractSteadyStateProblem,
    AbstractJumpProblem}

has_kwargs(_prob::AbstractDEProblem) = has_kwargs(typeof(_prob))
Base.@pure __has_kwargs(::Type{T}) where {T} = :kwargs ∈ fieldnames(T)
has_kwargs(::Type{T}) where {T} = __has_kwargs(T)

const allowedkeywords = (:dense,
    :saveat,
    :save_idxs,
    :tstops,
    :tspan,
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
    :progress_id,
    :timeseries_errors,
    :dense_errors,
    :weak_timeseries_errors,
    :weak_dense_errors,
    :wrap,
    :calculate_error,
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
    :u0,
    :p,
    # These two are from the default algorithm handling
    :default_set,
    :second_time,
    # This is for DiffEqDevTools
    :prob_choice,
    # Jump problems
    :alias_jump,
    # This is for copying/deepcopying noise in StochasticDiffEq
    :alias_noise,
    # This is for SimpleNonlinearSolve handling for batched Nonlinear Solves
    :batch,
    # Shooting method in BVP needs to differentiate between these two categories
    :nlsolve_kwargs,
    :odesolve_kwargs,
    # If Solvers which internally use linsolve
    :linsolve_kwargs,
    # Solvers internally using EnsembleProblem
    :ensemblealg)

const KWARGWARN_MESSAGE = """
                          Unrecognized keyword arguments found.
                          The only allowed keyword arguments to `solve` are:
                          $allowedkeywords

                          See https://diffeq.sciml.ai/stable/basics/common_solver_opts/ for more details.

                          Set kwargshandle=KeywordArgError for an error message.
                          Set kwargshandle=KeywordArgSilent to ignore this message.
                          """

const KWARGERROR_MESSAGE = """
                           Unrecognized keyword arguments found.
                           The only allowed keyword arguments to `solve` are:
                           $allowedkeywords

                           See https://diffeq.sciml.ai/stable/basics/common_solver_opts/ for more details.
                           """

struct CommonKwargError <: Exception
    kwargs::Any
end

function Base.showerror(io::IO, e::CommonKwargError)
    println(io, KWARGERROR_MESSAGE)
    notin = collect(map(x -> x ∉ allowedkeywords, keys(e.kwargs)))
    unrecognized = collect(keys(e.kwargs))[notin]
    print(io, "Unrecognized keyword arguments: ")
    printstyled(io, unrecognized; bold = true, color = :red)
    print(io, "\n\n")
    println(io, TruncatedStacktraces.VERBOSE_MSG)
end

@enum KeywordArgError KeywordArgWarn KeywordArgSilent

const INCOMPATIBLE_U0_MESSAGE = """
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
    println(io, TruncatedStacktraces.VERBOSE_MSG)
end

const NO_DEFAULT_ALGORITHM_MESSAGE = """
                                     Default algorithm choices require DifferentialEquations.jl.
                                     Please specify an algorithm (e.g., `solve(prob, Tsit5())` or
                                     `init(prob, Tsit5())` for an ODE) or import DifferentialEquations
                                     directly.

                                     You can find the list of available solvers at https://diffeq.sciml.ai/stable/solvers/ode_solve/
                                     and its associated pages.
                                     """

struct NoDefaultAlgorithmError <: Exception end

function Base.showerror(io::IO, e::NoDefaultAlgorithmError)
    print(io, NO_DEFAULT_ALGORITHM_MESSAGE)
    println(io, TruncatedStacktraces.VERBOSE_MSG)
end

const NO_TSPAN_MESSAGE = """
                         No tspan is set in the problem or chosen in the init/solve call
                         """

struct NoTspanError <: Exception end

function Base.showerror(io::IO, e::NoTspanError)
    print(io, NO_TSPAN_MESSAGE)
    println(io, TruncatedStacktraces.VERBOSE_MSG)
end

const NAN_TSPAN_MESSAGE = """
                          NaN tspan is set in the problem or chosen in the init/solve call.
                          Note that -Inf and Inf values are allowed in the timespan for solves
                          which are terminated via callbacks, however NaN values are not allowed
                          since the direction of time is undetermined.
                          """

struct NaNTspanError <: Exception end

function Base.showerror(io::IO, e::NaNTspanError)
    print(io, NAN_TSPAN_MESSAGE)
    println(io, TruncatedStacktraces.VERBOSE_MSG)
end

const NON_SOLVER_MESSAGE = """
                           The arguments to solve are incorrect.
                           The second argument must be a solver choice, `solve(prob,alg)`
                           where `alg` is a `<: AbstractDEAlgorithm`, e.g. `Tsit5()`.

                           Please double check the arguments being sent to the solver.

                           You can find the list of available solvers at https://diffeq.sciml.ai/stable/solvers/ode_solve/
                           and its associated pages.
                           """

struct NonSolverError <: Exception end

function Base.showerror(io::IO, e::NonSolverError)
    print(io, NON_SOLVER_MESSAGE)
    println(io, TruncatedStacktraces.VERBOSE_MSG)
end

const NOISE_SIZE_MESSAGE = """
                           Noise sizes are incompatible. The expected number of noise terms in the defined
                           `noise_rate_prototype` does not match the number of noise terms in the defined
                           `AbstractNoiseProcess`. Please ensure that
                           size(prob.noise_rate_prototype,2) == length(prob.noise.W[1]).

                           Note: Noise process definitions require that users specify `u0`, and this value is
                           directly used in the definition. For example, if `noise = WienerProcess(0.0,0.0)`,
                           then the noise process is a scalar with `u0=0.0`. If `noise = WienerProcess(0.0,[0.0])`,
                           then the noise process is a vector with `u0=0.0`. If `noise_rate_prototype = zeros(2,4)`,
                           then the noise process must be a 4-dimensional process, for example
                           `noise = WienerProcess(0.0,zeros(4))`. This error is a sign that the user definition
                           of `noise_rate_prototype` and `noise` are not aligned in this manner and the definitions should
                           be double checked.
                           """

struct NoiseSizeIncompatabilityError <: Exception
    prototypesize::Int
    noisesize::Int
end

function Base.showerror(io::IO, e::NoiseSizeIncompatabilityError)
    println(io, NOISE_SIZE_MESSAGE)
    println(io, "size(prob.noise_rate_prototype,2) = $(e.prototypesize)")
    println(io, "length(prob.noise.W[1]) = $(e.noisesize)")
    println(io, TruncatedStacktraces.VERBOSE_MSG)
end

const PROBSOLVER_PAIRING_MESSAGE = """
                                   Incompatible problem+solver pairing.
                                   For example, this can occur if an ODE solver is passed with an SDEProblem.
                                   Solvers are only capable of handling specific problem types. Please double
                                   check that the chosen pairing is capable for handling the given problems.
                                   """

struct ProblemSolverPairingError <: Exception
    prob::Any
    alg::Any
end

function Base.showerror(io::IO, e::ProblemSolverPairingError)
    println(io, PROBSOLVER_PAIRING_MESSAGE)
    println(io, "Problem type: $(SciMLBase.__parameterless_type(typeof(e.prob)))")
    println(io, "Solver type: $(SciMLBase.__parameterless_type(typeof(e.alg)))")
    println(io,
        "Problem types compatible with the chosen solver: $(compatible_problem_types(e.prob,e.alg))")
    println(io, TruncatedStacktraces.VERBOSE_MSG)
end

function compatible_problem_types(prob, alg)
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

const DIRECT_AUTODIFF_INCOMPATABILITY_MESSAGE = """
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
    println(io, TruncatedStacktraces.VERBOSE_MSG)
end

const NONCONCRETE_ELTYPE_MESSAGE = """
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
    eltype::Any
end

function Base.showerror(io::IO, e::NonConcreteEltypeError)
    print(io, NONCONCRETE_ELTYPE_MESSAGE)
    print(io, e.eltype)
    println(io, TruncatedStacktraces.VERBOSE_MSG)
end

const GENERIC_NUMBER_TYPE_ERROR_MESSAGE = """
                                          Non-standard number type (i.e. not Float32, Float64,
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
    alg::Any
    uType::Any
    tType::Any
end

function Base.showerror(io::IO, e::GenericNumberTypeError)
    println(io, GENERIC_NUMBER_TYPE_ERROR_MESSAGE)
    println(io, "Solver: $(e.alg)")
    println(io, "u0 type: $(e.uType)")
    print(io, "Timespan type: $(e.tType)")
    println(io, TruncatedStacktraces.VERBOSE_MSG)
end

const COMPLEX_SUPPORT_ERROR_MESSAGE = """
                                      Complex number type (i.e. ComplexF32, or ComplexF64)
                                      detected as the element type for the initial condition
                                      with an algorithm that does not support complex numbers.
                                      Please check that the initial condition type is correct.
                                      If complex number support is needed, try different solvers
                                      such as those from OrdinaryDiffEq.jl.
                                      """

struct ComplexSupportError <: Exception
    alg::Any
end

function Base.showerror(io::IO, e::ComplexSupportError)
    println(io, COMPLEX_SUPPORT_ERROR_MESSAGE)
    println(io, "Solver: $(e.alg)")
    println(io, TruncatedStacktraces.VERBOSE_MSG)
end

const COMPLEX_TSPAN_ERROR_MESSAGE = """
                                    Complex number type (i.e. ComplexF32, or ComplexF64)
                                    detected as the element type for the independent variable
                                    (i.e. time span). Please check that the tspan type is correct.
                                    No solvers support complex time spans. If this is required,
                                    please open an issue.
                                    """

struct ComplexTspanError <: Exception end

function Base.showerror(io::IO, e::ComplexTspanError)
    println(io, COMPLEX_TSPAN_ERROR_MESSAGE)
    println(io, TruncatedStacktraces.VERBOSE_MSG)
end

const TUPLE_STATE_ERROR_MESSAGE = """
                                  Tuple type used as a state. Since a tuple does not have vector
                                  properties, it will not work as a state type in equation solvers.
                                  Instead, change your equation from using tuple constructors `()`
                                  to static array constructors `SA[]`. For example, change:

                                  ```julia
                                  function ftup((a,b),p,t)
                                    return b,-a
                                  end
                                  u0 = (1.0,2.0)
                                  tspan = (0.0,1.0)
                                  ODEProblem(ftup,u0,tspan)
                                  ```

                                  to:

                                  ```julia
                                  using StaticArrays
                                  function fsa(u,p,t)
                                      SA[u[2],u[1]]
                                  end
                                  u0 = SA[1.0,2.0]
                                  tspan = (0.0,1.0)
                                  ODEProblem(ftup,u0,tspan)
                                  ```

                                  This will be safer and fast for small ODEs. For more information, see:
                                  https://diffeq.sciml.ai/stable/tutorials/faster_ode_example/#Further-Optimizations-of-Small-Non-Stiff-ODEs-with-StaticArrays
                                  """

struct TupleStateError <: Exception end

function Base.showerror(io::IO, e::TupleStateError)
    println(io, TUPLE_STATE_ERROR_MESSAGE)
    println(io, TruncatedStacktraces.VERBOSE_MSG)
end

const MASS_MATRIX_ERROR_MESSAGE = """
                                  Mass matrix size is incompatible with initial condition
                                  sizing. The mass matrix must represent the `vec`
                                  form of the initial condition `u0`, i.e.
                                  `size(mm,1) == size(mm,2) == length(u)`
                                  """

struct IncompatibleMassMatrixError <: Exception
    sz::Int
    len::Int
end

function Base.showerror(io::IO, e::IncompatibleMassMatrixError)
    println(io, MASS_MATRIX_ERROR_MESSAGE)
    print(io, "size(prob.f.mass_matrix,1): ")
    println(io, e.sz)
    print(io, "length(u0): ")
    println(e.len)
    println(io, TruncatedStacktraces.VERBOSE_MSG)
end

function init_call(_prob, args...; merge_callbacks = true, kwargshandle = nothing,
    kwargs...)
    kwargshandle = kwargshandle === nothing ? KeywordArgError : kwargshandle
    kwargshandle = has_kwargs(_prob) && haskey(_prob.kwargs, :kwargshandle) ?
                   _prob.kwargs[:kwargshandle] : kwargshandle

    if has_kwargs(_prob)
        if merge_callbacks && haskey(_prob.kwargs, :callback) && haskey(kwargs, :callback)
            kwargs_temp = NamedTuple{
                Base.diff_names(Base._nt_names(values(kwargs)),
                    (:callback,))}(values(kwargs))
            callbacks = NamedTuple{(:callback,)}((DiffEqBase.CallbackSet(_prob.kwargs[:callback],
                values(kwargs).callback),))
            kwargs = merge(kwargs_temp, callbacks)
        end
        kwargs = isempty(_prob.kwargs) ? kwargs : merge(values(_prob.kwargs), kwargs)
    end

    checkkwargs(kwargshandle; kwargs...)

    if _prob isa Union{ODEProblem, DAEProblem} && isnothing(_prob.u0)
        build_null_integrator(_prob, args...; kwargs...)
    elseif hasfield(typeof(_prob), :f) && hasfield(typeof(_prob.f), :f) &&
           typeof(_prob.f.f) <: EvalFunc
        Base.invokelatest(__init, _prob, args...; kwargs...)#::T
    else
        __init(_prob, args...; kwargs...)#::T
    end
end

function init(prob::Union{AbstractDEProblem, NonlinearProblem}, args...; sensealg = nothing,
    u0 = nothing, p = nothing, kwargs...)
    if sensealg === nothing && haskey(prob.kwargs, :sensealg)
        sensealg = prob.kwargs[:sensealg]
    end

    u0 = u0 !== nothing ? u0 : prob.u0
    p = p !== nothing ? p : prob.p

    init_up(prob, sensealg, u0, p, args...; kwargs...)
end

function init(prob::AbstractJumpProblem, args...; kwargs...)
    init_call(prob, args...; kwargs...)
end

function init_up(prob::AbstractDEProblem, sensealg, u0, p, args...; kwargs...)
    alg = extract_alg(args, kwargs, prob.kwargs)
    if isnothing(alg) # Default algorithm handling
        _prob = get_concrete_problem(prob, !(typeof(prob) <: DiscreteProblem); u0 = u0,
            p = p, kwargs...)
        init_call(_prob, args...; kwargs...)
    else
        _prob = get_concrete_problem(prob, isadaptive(alg); u0 = u0, p = p, kwargs...)
        _alg = prepare_alg(alg, _prob.u0, _prob.p, _prob)
        check_prob_alg_pairing(_prob, alg) # alg for improved inference
        if length(args) > 1
            init_call(_prob, _alg, Base.tail(args)...; kwargs...)
        else
            init_call(_prob, _alg; kwargs...)
        end
    end
end

function solve_call(_prob, args...; merge_callbacks = true, kwargshandle = nothing,
    kwargs...)
    kwargshandle = kwargshandle === nothing ? KeywordArgError : kwargshandle
    kwargshandle = has_kwargs(_prob) && haskey(_prob.kwargs, :kwargshandle) ?
                   _prob.kwargs[:kwargshandle] : kwargshandle

    if has_kwargs(_prob)
        if merge_callbacks && haskey(_prob.kwargs, :callback) && haskey(kwargs, :callback)
            kwargs_temp = NamedTuple{
                Base.diff_names(Base._nt_names(values(kwargs)),
                    (:callback,))}(values(kwargs))
            callbacks = NamedTuple{(:callback,)}((DiffEqBase.CallbackSet(_prob.kwargs[:callback],
                values(kwargs).callback),))
            kwargs = merge(kwargs_temp, callbacks)
        end
        kwargs = isempty(_prob.kwargs) ? kwargs : merge(values(_prob.kwargs), kwargs)
    end

    checkkwargs(kwargshandle; kwargs...)
    if isdefined(_prob, :u0)
        if _prob.u0 isa Array &&
           !isconcretetype(RecursiveArrayTools.recursive_unitless_eltype(_prob.u0))
            throw(NonConcreteEltypeError(RecursiveArrayTools.recursive_unitless_eltype(_prob.u0)))
        end

        if _prob.u0 === nothing
            return build_null_solution(_prob, args...; kwargs...)
        end
    end

    if hasfield(typeof(_prob), :f) && hasfield(typeof(_prob.f), :f) &&
       typeof(_prob.f.f) <: EvalFunc
        Base.invokelatest(__solve, _prob, args...; kwargs...)#::T
    else
        __solve(_prob, args...; kwargs...)#::T
    end
end

mutable struct NullODEIntegrator{IIP, ProbType, T, SolType, F, P} <:
               AbstractODEIntegrator{Nothing, IIP, Nothing, T}
    du::Vector{Float64}
    u::Vector{Float64}
    t::T
    prob::ProbType
    sol::SolType
    f::F
    p::P
end
function build_null_integrator(prob::AbstractDEProblem, args...;
    kwargs...)
    sol = solve(prob, args...; kwargs...)
    return NullODEIntegrator{isinplace(prob), typeof(prob), eltype(prob.tspan), typeof(sol),
        typeof(prob.f), typeof(prob.p),
    }(Float64[],
        Float64[],
        first(prob.tspan),
        prob,
        sol,
        prob.f,
        prob.p)
end
function solve!(integ::NullODEIntegrator)
    integ.t = integ.sol.t[end]
    return nothing
end
function step!(integ::NullODEIntegrator, dt = nothing, stop_at_tdt = false)
    if !isnothing(dt)
        integ.t += dt
    else
        integ.t = integ.sol[end]
    end
    return nothing
end

function build_null_solution(prob::AbstractDEProblem, args...;
    saveat = (),
    save_everystep = true,
    save_on = true,
    save_start = save_everystep || isempty(saveat) ||
                     saveat isa Number || prob.tspan[1] in saveat,
    save_end = true,
    kwargs...)
    ts = if saveat === ()
        if save_start && save_end
            [prob.tspan[1], prob.tspan[2]]
        elseif save_start && !save_end
            [prob.tspan[1]]
        elseif !save_start && save_end
            [prob.tspan[2]]
        else
            eltype(prob.tspan)[]
        end
    elseif saveat isa Number
        prob.tspan[1]:saveat:prob.tspan[2]
    else
        saveat
    end

    timeseries = [Float64[] for i in 1:length(ts)]

    build_solution(prob, nothing, ts, timeseries, retcode = ReturnCode.Success)
end

function build_null_solution(prob::Union{SteadyStateProblem, NonlinearProblem}, args...;
    saveat = (),
    save_everystep = true,
    save_on = true,
    save_start = save_everystep || isempty(saveat) ||
                     saveat isa Number || prob.tspan[1] in saveat,
    save_end = true,
    kwargs...)
    SciMLBase.build_solution(prob, nothing, Float64[], nothing;
        retcode = ReturnCode.Success)
end

"""
```julia
solve(prob::AbstractDEProblem, alg::Union{AbstractDEAlgorithm,Nothing}; kwargs...)
```

## Arguments

The only positional argument is `alg` which is optional. By default, `alg = nothing`.
If `alg = nothing`, then `solve` dispatches to the DifferentialEquations.jl automated
algorithm selection (if `using DifferentialEquations` was done, otherwise it will
error with a `MethodError`).

## Keyword Arguments

The DifferentialEquations.jl universe has a large set of common arguments available
for the `solve` function. These arguments apply to `solve` on any problem type and
are only limited by limitations of the specific implementations.

Many of the defaults depend on the algorithm or the package the algorithm derives
from. Not all of the interface is provided by every algorithm.
For more detailed information on the defaults and the available options
for specific algorithms / packages, see the manual pages for the solvers of specific
problems. To see whether a specific package is compatible with the use of a
given option, see the [Solver Compatibility Chart](https://docs.sciml.ai/DiffEqDocs/stable/basics/compatibility_chart/#Solver-Compatibility-Chart)

### Default Algorithm Hinting

To help choose the default algorithm, the keyword argument `alg_hints` is
provided to `solve`. `alg_hints` is a `Vector{Symbol}` which describe the
problem at a high level to the solver. The options are:

* `:auto` vs `:nonstiff` vs `:stiff` - Denotes the equation as nonstiff/stiff.
  `:auto` allow the default handling algorithm to choose stiffness detection
  algorithms. The default handling defaults to using `:auto`.

Currently unused options include:

* `:interpolant` - Denotes that a high-precision interpolation is important.
* `:memorybound` - Denotes that the solver will be memory bound.

This functionality is derived via the benchmarks in
[SciMLBenchmarks.jl](https://github.com/SciML/SciMLBenchmarks.jl)

#### SDE Specific Alghints

* `:additive` - Denotes that the underlying SDE has additive noise.
* `:stratonovich` - Denotes that the solution should adhere to the Stratonovich
  interpretation.

### Output Control

These arguments control the output behavior of the solvers. It defaults to maximum
output to give the best interactive user experience, but can be reduced all the
way to only saving the solution at the final timepoint.

The following options are all related to output control. See the "Examples"
section at the end of this page for some example usage.

* `dense`: Denotes whether to save the extra pieces required for dense (continuous)
  output. Default is `save_everystep && isempty(saveat)` for algorithms which have
  the ability to produce dense output, i.e. by default it's `true` unless the user
  has turned off saving on steps or has chosen a `saveat` value. If `dense=false`,
  the solution still acts like a function, and `sol(t)` is a linear interpolation
  between the saved time points.
* `saveat`: Denotes specific times to save the solution at, during the solving
  phase. The solver will save at each of the timepoints in this array in the
  most efficient manner available to the solver. If only `saveat` is given, then
  the arguments `save_everystep` and `dense` are `false` by default.
  If `saveat` is given a number, then it will automatically expand to
  `tspan[1]:saveat:tspan[2]`. For methods where interpolation is not possible,
  `saveat` may be equivalent to `tstops`. The default value is `[]`.
* `save_idxs`: Denotes the indices for the components of the equation to save.
  Defaults to saving all indices. For example, if you are solving a 3-dimensional ODE,
  and given `save_idxs = [1, 3]`, only the first and third components of the
  solution will be outputted.
  Notice that of course in this case the outputted solution will be two-dimensional.
* `tstops`: Denotes *extra* times that the timestepping algorithm must step to.
  This should be used to help the solver deal with discontinuities and
  singularities, since stepping exactly at the time of the discontinuity will
  improve accuracy. If a method cannot change timesteps (fixed timestep
  multistep methods), then `tstops` will use an interpolation,
  matching the behavior of `saveat`. If a method cannot change timesteps and
  also cannot interpolate, then `tstops` must be a multiple of `dt` or else an
  error will be thrown. Default is `[]`.
* `d_discontinuities:` Denotes locations of discontinuities in low order derivatives.
  This will force FSAL algorithms which assume derivative continuity to re-evaluate
  the derivatives at the point of discontinuity. The default is `[]`.
* `save_everystep`: Saves the result at every step.
  Default is true if `isempty(saveat)`.
* `save_on`: Denotes whether intermediate solutions are saved. This overrides the
  settings of `dense`, `saveat` and `save_everystep` and is used by some applications
  to manually turn off saving temporarily. Everyday use of the solvers should leave
  this unchanged. Defaults to `true`.
* `save_start`: Denotes whether the initial condition should be included in
  the solution type as the first timepoint. Defaults to `true`.
* `save_end`: Denotes whether the final timepoint is forced to be saved,
  regardless of the other saving settings. Defaults to `true`.
* `initialize_save`: Denotes whether to save after the callback initialization
  phase (when `u_modified=true`). Defaults to `true`.

Note that `dense` requires `save_everystep=true` and `saveat=false`. If you need
additional saving while keeping dense output, see
[the SavingCallback in the Callback Library](https://docs.sciml.ai/DiffEqCallbacks/stable/output_saving/#DiffEqCallbacks.SavingCallback).

### Stepsize Control

These arguments control the timestepping routines.

#### Basic Stepsize Control

These are the standard options for controlling stepping behavior. Error estimates
do the comparison

```math
err_{scaled} = err/(abstol + max(uprev,u)*reltol)
```

The scaled error is guaranteed to be `<1` for a given local error estimate
(note: error estimates are local unless the method specifies otherwise). `abstol`
controls the non-scaling error and thus can be thought of as the error around zero.
`reltol` scales with the size of the dependent variables and so one can interpret
`reltol=1e-3` as roughly being (locally) correct to 3 digits. Note tolerances can
be specified element-wise by passing a vector whose size matches `u0`.

* `adaptive`: Turns on adaptive timestepping for appropriate methods. Default
  is true.
* `abstol`: Absolute tolerance in adaptive timestepping. This is the tolerance
  on local error estimates, not necessarily the global error (though these quantities
  are related). Defaults to `1e-6` on deterministic equations (ODEs/DDEs/DAEs) and `1e-2`
  on stochastic equations (SDEs/RODEs).
* `reltol`: Relative tolerance in adaptive timestepping.  This is the tolerance
  on local error estimates, not necessarily the global error (though these quantities
  are related). Defaults to `1e-3` on deterministic equations (ODEs/DDEs/DAEs) and `1e-2`
  on stochastic equations (SDEs/RODEs).
* `dt`: Sets the initial stepsize. This is also the stepsize for fixed
  timestep methods. Defaults to an automatic choice if the method is adaptive.
* `dtmax`: Maximum dt for adaptive timestepping. Defaults are
  package-dependent.
* `dtmin`: Minimum dt for adaptive timestepping. Defaults are
  package-dependent.
* `force_dtmin`: Declares whether to continue, forcing the minimum `dt` usage.
  Default is `false`, which has the solver throw a warning and exit early when
  encountering the minimum `dt`. Setting this true allows the solver to continue,
  never letting `dt` go below `dtmin` (and ignoring error tolerances in those
  cases). Note that `true` is not compatible with most interop packages.

#### Fixed Stepsize Usage

Note that if a method does not have adaptivity, the following rules apply:

* If `dt` is set, then the algorithm will step with size `dt` each iteration.
* If `tstops` and `dt` are both set, then the algorithm will step with either a
  size `dt`, or use a smaller step to hit the `tstops` point.
* If `tstops` is set without `dt`, then the algorithm will step directly to
  each value in `tstops`
* If neither `dt` nor `tstops` are set, the solver will throw an error.

#### [Advanced Adaptive Stepsize Control](https://docs.sciml.ai/DiffEqDocs/stable/extras/timestepping/)

These arguments control more advanced parts of the internals of adaptive timestepping
and are mostly used to make it more efficient on specific problems. For detained
explanations of the timestepping algorithms, see the
[timestepping descriptions](https://docs.sciml.ai/DiffEqDocs/stable/extras/timestepping/#timestepping)

* `internalnorm`: The norm function `internalnorm(u,t)` which error estimates
  are calculated. Required are two dispatches: one dispatch for the state variable
  and the other on the elements of the state variable (scalar norm).
  Defaults are package-dependent.
* `controller`: Possible examples are [`IController`](https://docs.sciml.ai/DiffEqDocs/stable/extras/timestepping/#OrdinaryDiffEq.IController),
  [`PIController`](https://docs.sciml.ai/DiffEqDocs/stable/extras/timestepping/#OrdinaryDiffEq.PIController),
  [`PIDController`](https://docs.sciml.ai/DiffEqDocs/stable/extras/timestepping/#OrdinaryDiffEq.PIDController),
  [`PredictiveController`](https://docs.sciml.ai/DiffEqDocs/stable/extras/timestepping/#OrdinaryDiffEq.PredictiveController).
  Default is algorithm-dependent.
* `gamma`: The risk-factor γ in the q equation for adaptive timestepping
  of the controllers using it.
  Default is algorithm-dependent.
* `beta1`: The Lund stabilization α parameter.
  Default is algorithm-dependent.
* `beta2`: The Lund stabilization β parameter.
  Default is algorithm-dependent.
* `qmax`: Defines the maximum value possible for the adaptive q.
  Default is algorithm-dependent.
* `qmin`: Defines the minimum value possible for the adaptive q.
  Default is algorithm-dependent.
* `qsteady_min`: Defines the minimum for the range around 1 where the timestep
  is held constant. Default is algorithm-dependent.
* `qsteady_max`: Defines the maximum for the range around 1 where the timestep
  is held constant. Default is algorithm-dependent.
* `qoldinit`: The initial `qold` in stabilization stepping.
  Default is algorithm-dependent.
* `failfactor`: The amount to decrease the timestep by if the Newton iterations
  of an implicit method fail. Default is 2.

### Memory Optimizations

* `calck`: Turns on and off the internal ability for intermediate
  interpolations (also known as intermediate density). Not the same as `dense`, which is post-solution interpolation.
  This defaults to `dense || !isempty(saveat) ||  "no custom callback is given"`.
  This can be used to turn off interpolations
  (to save memory) if one isn't using interpolations when a custom callback is
  used. Another case where this may be used is to turn on interpolations for
  usage in the integrator interface even when interpolations are used nowhere else.
  Note that this is only required if the algorithm doesn't have
  a free or lazy interpolation (`DP8()`). If `calck = false`, `saveat` cannot be used.
  The rare keyword `calck` can be useful in event handling.
* `alias_u0`: allows the solver to alias the initial condition array that is contained
  in the problem struct. Defaults to false.

### Miscellaneous

* `maxiters`: Maximum number of iterations before stopping. Defaults to 1e5.
* `callback`: Specifies a callback. Defaults to a callback function which
  performs the saving routine. For more information, see the
  [Event Handling and Callback Functions manual page](https://docs.sciml.ai/DiffEqCallbacks/stable/).
* `isoutofdomain`: Specifies a function `isoutofdomain(u,p,t)` where, when it
  returns true, it will reject the timestep. Disabled by default.
* `unstable_check`: Specifies a function `unstable_check(dt,u,p,t)` where, when
  it returns true, it will cause the solver to exit and throw a warning. Defaults
  to `any(isnan,u)`, i.e. checking if any value is a NaN.
* `verbose`: Toggles whether warnings are thrown when the solver exits early.
  Defaults to true.
* `merge_callbacks`: Toggles whether to merge `prob.callback` with the `solve` keyword
  argument `callback`. Defaults to `true`.
* `wrap`: Toggles whether to wrap the solution if `prob.problem_type` has a preferred
  alternate wrapper type for the solution. Useful when speed, but not shape of solution
  is important. Defaults to `Val(true)`. `Val(false)` will cancel wrapping the solution.

### Progress Monitoring

These arguments control the usage of the progressbar in ProgressLogging.jl compatible environments.

* `progress`: Turns on/off the Juno progressbar. Default is false.
* `progress_steps`: Numbers of steps between updates of the progress bar.
  Default is 1000.
* `progress_name`: Controls the name of the progressbar. Default is the name
  of the problem type.
* `progress_message`: Controls the message with the progressbar. Defaults to
  showing `dt`, `t`, the maximum of `u`.
* `progress_id`: Controls the ID of the progress log message to distinguish simultaneous simulations.

### Error Calculations

If you are using the test problems (ex: `ODETestProblem`), then the following
options control the errors which are calculated:

* `timeseries_errors`: Turns on and off the calculation of errors at the steps
  which were taken, such as the `l2` error. Default is true.
* `dense_errors`: Turns on and off the calculation of errors at the steps which
  require dense output and calculate the error at 100 evenly-spaced points
  throughout `tspan`. An example is the `L2` error. Default is false.

### Sensitivity Algorithms (`sensealg`)

`sensealg` is used for choosing the way the automatic differentiation is performed.
For more information, see the documentation for SciMLSensitivity:
https://docs.sciml.ai/SciMLSensitivity/stable/

## Examples

The following lines are examples of how one could use the configuration of
`solve()`. For these examples a 3-dimensional ODE problem is assumed, however
the extension to other types is straightforward.

1. `solve(prob, AlgorithmName())` : The "default" setting, with a user-specified
  algorithm (given by `AlgorithmName()`). All parameters get their default values.
  This means that the solution is saved at the steps the Algorithm stops internally
  and dense output is enabled if the chosen algorithm allows for it.

  All other integration parameters (e.g. stepsize) are chosen automatically.
2. `solve(prob, saveat = 0.01, abstol = 1e-9, reltol = 1e-9)` : Standard setting
  for accurate output at specified (and equidistant) time intervals, used for
  e.g. Fourier Transform. The solution is given every 0.01 time units,
  starting from `tspan[1]`. The solver used is `Tsit5()` since no keyword
  `alg_hits` is given.

3. `solve(prob, maxiters = 1e7, progress = true, save_idxs = [1])` : Using longer
  maximum number of solver iterations can be useful when a given `tspan` is very
  long. This example only saves the first of the variables of the system, either
  to save size or because the user does not care about the others. Finally, with
  `progress = true` you are enabling the progress bar.
"""
function solve(prob::AbstractDEProblem, args...; sensealg = nothing,
    u0 = nothing, p = nothing, wrap = Val(true), kwargs...)
    if sensealg === nothing && haskey(prob.kwargs, :sensealg)
        sensealg = prob.kwargs[:sensealg]
    end

    u0 = u0 !== nothing ? u0 : prob.u0
    p = p !== nothing ? p : prob.p

    if wrap isa Val{true}
        wrap_sol(solve_up(prob, sensealg, u0, p, args...; kwargs...))
    else
        solve_up(prob, sensealg, u0, p, args...; kwargs...)
    end
end

"""
```julia
solve(prob::NonlinearProblem, alg::Union{AbstractNonlinearAlgorithm,Nothing}; kwargs...)
```

## Arguments

The only positional argument is `alg` which is optional. By default, `alg = nothing`.
If `alg = nothing`, then `solve` dispatches to the NonlinearSolve.jl automated
algorithm selection (if `using NonlinearSolve` was done, otherwise it will
error with a `MethodError`).

## Keyword Arguments

The NonlinearSolve.jl universe has a large set of common arguments available
for the `solve` function. These arguments apply to `solve` on any problem type and
are only limited by limitations of the specific implementations.

Many of the defaults depend on the algorithm or the package the algorithm derives
from. Not all of the interface is provided by every algorithm.
For more detailed information on the defaults and the available options
for specific algorithms / packages, see the manual pages for the solvers of specific
problems.

#### Error Control

* `abstol`: Absolute tolerance.
* `reltol`: Relative tolerance.

### Miscellaneous

* `maxiters`: Maximum number of iterations before stopping. Defaults to 1e5.
* `verbose`: Toggles whether warnings are thrown when the solver exits early.
  Defaults to true.

### Sensitivity Algorithms (`sensealg`)

`sensealg` is used for choosing the way the automatic differentiation is performed.
    For more information, see the documentation for SciMLSensitivity:
    https://docs.sciml.ai/SciMLSensitivity/stable/
"""
function solve(prob::NonlinearProblem, args...; sensealg = nothing,
    u0 = nothing, p = nothing, wrap = Val(true), kwargs...)
    if sensealg === nothing && haskey(prob.kwargs, :sensealg)
        sensealg = prob.kwargs[:sensealg]
    end

    u0 = u0 !== nothing ? u0 : prob.u0
    p = p !== nothing ? p : prob.p

    if wrap isa Val{true}
        wrap_sol(solve_up(prob, sensealg, u0, p, args...; kwargs...))
    else
        solve_up(prob, sensealg, u0, p, args...; kwargs...)
    end
end

function solve_up(prob::Union{AbstractDEProblem, NonlinearProblem}, sensealg, u0, p,
    args...; kwargs...)
    alg = extract_alg(args, kwargs, prob.kwargs)
    if isnothing(alg) # Default algorithm handling
        _prob = get_concrete_problem(prob, !(typeof(prob) <: DiscreteProblem); u0 = u0,
            p = p, kwargs...)
        solve_call(_prob, args...; kwargs...)
    else
        _prob = get_concrete_problem(prob, isadaptive(alg); u0 = u0, p = p, kwargs...)
        _alg = prepare_alg(alg, _prob.u0, _prob.p, _prob)
        check_prob_alg_pairing(_prob, alg) # use alg for improved inference
        if length(args) > 1
            solve_call(_prob, _alg, Base.tail(args)...; kwargs...)
        else
            solve_call(_prob, _alg; kwargs...)
        end
    end
end

function solve_call(prob::SteadyStateProblem,
    alg::SciMLBase.AbstractNonlinearAlgorithm, args...;
    kwargs...)
    solve_call(NonlinearProblem(prob),
        alg, args...;
        kwargs...)
end

function solve(prob::EnsembleProblem, args...; kwargs...)
    if isempty(args) || length(args) == 1 && typeof(args[1]) <: EnsembleAlgorithm
        __solve(prob, nothing, args...; kwargs...)
    else
        __solve(prob, args...; kwargs...)
    end
end
function solve(prob::SciMLBase.WeightedEnsembleProblem, args...; kwargs...)
    SciMLBase.WeightedEnsembleSolution(solve(prob.ensembleprob), prob.weights)
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
            unrecognized = setdiff(keys(kwargs), allowedkeywords)
            print("Unrecognized keyword arguments: ")
            printstyled(unrecognized; bold = true, color = :red)
            print("\n\n")
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

function solve(prob::PDEProblem, alg::AbstractDEAlgorithm, args...;
    kwargs...)
    solve(prob.prob, alg, args...; kwargs...)
end

function init(prob::PDEProblem, alg::AbstractDEAlgorithm, args...;
    kwargs...)
    init(prob.prob, alg, args...; kwargs...)
end

function get_concrete_problem(prob, isadapt; kwargs...)
    p = get_concrete_p(prob, kwargs)
    tspan = get_concrete_tspan(prob, isadapt, kwargs, p)
    u0 = get_concrete_u0(prob, isadapt, tspan[1], kwargs)
    u0_promote = promote_u0(u0, p, tspan[1])
    tspan_promote = promote_tspan(u0_promote, p, tspan, prob, kwargs)
    f_promote = promote_f(prob.f, Val(SciMLBase.specialization(prob.f)), u0_promote, p,
        tspan_promote[1])
    if isconcreteu0(prob, tspan[1], kwargs) && typeof(u0_promote) === typeof(prob.u0) &&
       prob.tspan == tspan && typeof(prob.tspan) === typeof(tspan_promote) &&
       p === prob.p && f_promote === prob.f
        return prob
    else
        return remake(prob; f = f_promote, u0 = u0_promote, p = p, tspan = tspan_promote)
    end
end

function get_concrete_problem(prob::DAEProblem, isadapt; kwargs...)
    p = get_concrete_p(prob, kwargs)
    tspan = get_concrete_tspan(prob, isadapt, kwargs, p)
    u0 = get_concrete_u0(prob, isadapt, tspan[1], kwargs)
    du0 = get_concrete_du0(prob, isadapt, tspan[1], kwargs)

    u0_promote = promote_u0(u0, p, tspan[1])
    du0_promote = promote_u0(du0, p, tspan[1])
    tspan_promote = promote_tspan(u0_promote, p, tspan, prob, kwargs)

    f_promote = promote_f(prob.f, Val(SciMLBase.specialization(prob.f)), u0_promote, p,
        tspan_promote[1])
    if isconcreteu0(prob, tspan[1], kwargs) && typeof(u0_promote) === typeof(prob.u0) &&
       isconcretedu0(prob, tspan[1], kwargs) && typeof(du0_promote) === typeof(prob.du0) &&
       prob.tspan == tspan && typeof(prob.tspan) === typeof(tspan_promote) &&
       p === prob.p && f_promote === prob.f
        return prob
    else
        return remake(prob; f = f_promote, du0 = du0_promote, u0 = u0_promote, p = p,
            tspan = tspan_promote)
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

    remake(prob; u0 = u0, tspan = tspan, p = p, constant_lags = constant_lags)
end

function promote_f(f::F, ::Val{specialize}, u0, p, t) where {F, specialize}
    # Ensure our jacobian will be of the same type as u0
    uElType = u0 === nothing ? Float64 : eltype(u0)
    if isdefined(f, :jac_prototype) && f.jac_prototype isa AbstractArray
        f = @set f.jac_prototype = similar(f.jac_prototype, uElType)
    end

    @static if VERSION >= v"1.8-"
        f = if f isa ODEFunction && isinplace(f) && !(f.f isa AbstractSciMLOperator) &&
               # Some reinitialization code still uses NLSolvers stuff which doesn't
               # properly tag, so opt-out if potentially a mass matrix DAE
               f.mass_matrix isa UniformScaling &&
               # Jacobians don't wrap, so just ignore those cases
               f.jac === nothing &&
               ((specialize === SciMLBase.AutoSpecialize && eltype(u0) !== Any &&
                 RecursiveArrayTools.recursive_unitless_eltype(u0) === eltype(u0) &&
                 one(t) === oneunit(t) &&
                 Tricks.static_hasmethod(ArrayInterface.promote_eltype,
                     Tuple{Type{typeof(u0)}, Type{dualgen(eltype(u0))}}) &&
                 Tricks.static_hasmethod(promote_rule,
                     Tuple{Type{eltype(u0)}, Type{dualgen(eltype(u0))}}) &&
                 Tricks.static_hasmethod(promote_rule,
                     Tuple{Type{eltype(u0)}, Type{typeof(t)}})) ||
                (specialize === SciMLBase.FunctionWrapperSpecialize &&
                 !(f.f isa FunctionWrappersWrappers.FunctionWrappersWrapper)))
            return unwrapped_f(f, wrapfun_iip(f.f, (u0, u0, p, t)))
        else
            return f
        end
    else
        return f
    end
end

function promote_f(f::SplitFunction, ::Val{specialize}, u0, p, t) where {specialize}
    typeof(f.cache) === typeof(u0) && isinplace(f) ? f : remake(f, cache = zero(u0))
end
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

    any(isnan, tspan) && throw(NaNTspanError())

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

    nu0 = length(something(_u0, ()))
    if isdefined(prob.f, :mass_matrix) && prob.f.mass_matrix !== nothing &&
       prob.f.mass_matrix isa AbstractArray &&
       size(prob.f.mass_matrix, 1) !== nu0
        throw(IncompatibleMassMatrixError(size(prob.f.mass_matrix, 1), nu0))
    end

    if _u0 isa Tuple
        throw(TupleStateError())
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

eval_u0(u0::Function) = true
eval_u0(u0) = false

function __solve(prob::AbstractDEProblem, args...; default_set = false, second_time = false,
    kwargs...)
    if second_time
        throw(NoDefaultAlgorithmError())
    elseif length(args) > 0 && !(typeof(args[1]) <: Union{Nothing, AbstractDEAlgorithm})
        throw(NonSolverError())
    else
        __solve(prob, nothing, args...; default_set = false, second_time = true, kwargs...)
    end
end

function __init(prob::AbstractDEProblem, args...; default_set = false, second_time = false,
    kwargs...)
    if second_time
        throw(NoDefaultAlgorithmError())
    elseif length(args) > 0 && !(typeof(args[1]) <: Union{Nothing, AbstractDEAlgorithm})
        throw(NonSolverError())
    else
        __init(prob, nothing, args...; default_set = false, second_time = true, kwargs...)
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

    if isdefined(prob, :u0) && eltype(prob.u0) <: ForwardDiff.Dual &&
       !SciMLBase.isautodifferentiable(alg)
        throw(DirectAutodiffError())
    end

    if prob isa SDEProblem && prob.noise_rate_prototype !== nothing &&
       prob.noise !== nothing &&
       size(prob.noise_rate_prototype, 2) != length(prob.noise.W[1])
        throw(NoiseSizeIncompatabilityError(size(prob.noise_rate_prototype, 2),
            length(prob.noise.W[1])))
    end

    # Complex number support comes before arbitrary number support for a more direct
    # error message.
    if !SciMLBase.allowscomplex(alg)
        if isdefined(prob, :u0) &&
           RecursiveArrayTools.recursive_unitless_eltype(prob.u0) <: Complex
            throw(ComplexSupportError(alg))
        end
    end

    if isdefined(prob, :tspan) && eltype(prob.tspan) <: Complex
        throw(ComplexTspanError())
    end

    # Check for concrete element type so that the non-concrete case throws a better error
    if !SciMLBase.allows_arbitrary_number_types(alg)
        if isdefined(prob, :u0)
            uType = RecursiveArrayTools.recursive_unitless_eltype(prob.u0)
            if Base.isconcretetype(uType) &&
               !(uType <: Union{Float32, Float64, ComplexF32, ComplexF64})
                throw(GenericNumberTypeError(alg,
                    isdefined(prob, :u0) ? typeof(prob.u0) :
                    nothing,
                    isdefined(prob, :tspan) ? typeof(prob.tspan) :
                    nothing))
            end
        end

        if isdefined(prob, :tspan)
            tType = eltype(prob.tspan)
            if Base.isconcretetype(tType) &&
               !(tType <: Union{Float32, Float64, ComplexF32, ComplexF64})
                throw(GenericNumberTypeError(alg,
                    isdefined(prob, :u0) ? typeof(prob.u0) :
                    nothing,
                    isdefined(prob, :tspan) ? typeof(prob.tspan) :
                    nothing))
            end
        end
    end
end

################### Differentiation

"""
Ignores all adjoint definitions (i.e. `sensealg`) and proceeds to do standard
AD through the `solve` functions. Generally only used internally for implementing
discrete sensitivity algorithms.
"""
struct SensitivityADPassThrough <: AbstractDEAlgorithm end

function ChainRulesCore.frule(::typeof(solve_up), prob,
    sensealg::Union{Nothing, AbstractSensitivityAlgorithm},
    u0, p, args...;
    kwargs...)
    _solve_forward(prob, sensealg, u0, p, SciMLBase.ChainRulesOriginator(), args...;
        kwargs...)
end

function ChainRulesCore.rrule(::typeof(solve_up), prob::AbstractDEProblem,
    sensealg::Union{Nothing, AbstractSensitivityAlgorithm},
    u0, p, args...;
    kwargs...)
    _solve_adjoint(prob, sensealg, u0, p, SciMLBase.ChainRulesOriginator(), args...;
        kwargs...)
end

###
### Legacy Dispatches to be Non-Breaking
###

@deprecate concrete_solve(prob::AbstractDEProblem,
    alg::Union{AbstractDEAlgorithm, Nothing},
    u0 = prob.u0, p = prob.p, args...; kwargs...) solve(prob, alg,
    args...;
    u0 = u0,
    p = p,
    kwargs...)

function _solve_adjoint(prob, sensealg, u0, p, originator, args...; merge_callbacks = true,
    kwargs...)
    alg = extract_alg(args, kwargs, prob.kwargs)
    if isnothing(alg) # Default algorithm handling
        _prob = get_concrete_problem(prob, !(typeof(prob) <: DiscreteProblem); u0 = u0,
            p = p, kwargs...)
    else
        _prob = get_concrete_problem(prob, isadaptive(alg); u0 = u0, p = p, kwargs...)
    end

    if has_kwargs(_prob)
        if merge_callbacks && haskey(_prob.kwargs, :callback) && haskey(kwargs, :callback)
            kwargs_temp = NamedTuple{
                Base.diff_names(Base._nt_names(values(kwargs)),
                    (:callback,))}(values(kwargs))
            callbacks = NamedTuple{(:callback,)}((DiffEqBase.CallbackSet(_prob.kwargs[:callback],
                values(kwargs).callback),))
            kwargs = merge(kwargs_temp, callbacks)
        end
        kwargs = isempty(_prob.kwargs) ? kwargs : merge(values(_prob.kwargs), kwargs)
    end

    if isempty(args)
        _concrete_solve_adjoint(_prob, alg, sensealg, u0, p, originator; kwargs...)
    else
        _concrete_solve_adjoint(_prob, alg, sensealg, u0, p, originator,
            Base.tail(args)...; kwargs...)
    end
end

function _solve_forward(prob, sensealg, u0, p, originator, args...; merge_callbacks = true,
    kwargs...)
    alg = extract_alg(args, kwargs, prob.kwargs)
    if isnothing(alg) # Default algorithm handling
        _prob = get_concrete_problem(prob, !(typeof(prob) <: DiscreteProblem); u0 = u0,
            p = p, kwargs...)
    else
        _prob = get_concrete_problem(prob, isadaptive(alg); u0 = u0, p = p, kwargs...)
    end

    if has_kwargs(_prob)
        if merge_callbacks && haskey(_prob.kwargs, :callback) && haskey(kwargs, :callback)
            kwargs_temp = NamedTuple{
                Base.diff_names(Base._nt_names(values(kwargs)),
                    (:callback,))}(values(kwargs))
            callbacks = NamedTuple{(:callback,)}((DiffEqBase.CallbackSet(_prob.kwargs[:callback],
                values(kwargs).callback),))
            kwargs = merge(kwargs_temp, callbacks)
        end
        kwargs = isempty(_prob.kwargs) ? kwargs : merge(values(_prob.kwargs), kwargs)
    end

    if isempty(args)
        _concrete_solve_forward(_prob, alg, sensealg, u0, p, originator; kwargs...)
    else
        _concrete_solve_forward(_prob, alg, sensealg, u0, p, originator,
            Base.tail(args)...; kwargs...)
    end
end

@inline function extract_alg(solve_args, solve_kwargs, prob_kwargs)
    if isempty(solve_args) || isnothing(solve_args[1])
        if haskey(solve_kwargs, :alg)
            solve_kwargs[:alg]
        elseif haskey(prob_kwargs, :alg)
            prob_kwargs[:alg]
        else
            nothing
        end
    elseif solve_args[1] isa SciMLBase.AbstractSciMLAlgorithm
        solve_args[1]
    else
        nothing
    end
end

####
# Catch undefined AD overload cases

const ADJOINT_NOT_FOUND_MESSAGE = """
                                  Compatibility with reverse-mode automatic differentiation requires SciMLSensitivity.jl.
                                  Please install SciMLSensitivity.jl and do `using SciMLSensitivity`/`import SciMLSensitivity`
                                  for this functionality. For more details, see https://sensitivity.sciml.ai/dev/.
                                  """

struct AdjointNotFoundError <: Exception end

function Base.showerror(io::IO, e::AdjointNotFoundError)
    print(io, ADJOINT_NOT_FOUND_MESSAGE)
    println(io, TruncatedStacktraces.VERBOSE_MSG)
end

function _concrete_solve_adjoint(args...; kwargs...)
    throw(AdjointNotFoundError())
end

const FORWARD_SENSITIVITY_NOT_FOUND_MESSAGE = """
                                              Compatibility with forward-mode automatic differentiation requires SciMLSensitivity.jl.
                                              Please install SciMLSensitivity.jl and do `using SciMLSensitivity`/`import SciMLSensitivity`
                                              for this functionality. For more details, see https://sensitivity.sciml.ai/dev/.
                                              """

struct ForwardSensitivityNotFoundError <: Exception end

function Base.showerror(io::IO, e::ForwardSensitivityNotFoundError)
    print(io, FORWARD_SENSITIVITY_NOT_FOUND_MESSAGE)
    println(io, TruncatedStacktraces.VERBOSE_MSG)
end

function _concrete_solve_forward(args...; kwargs...)
    throw(ForwardSensitivityNotFoundError())
end
