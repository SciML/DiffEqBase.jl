"""
    NonlinearSafeTerminationReturnCode

Return Codes for the safe nonlinear termination conditions.
"""
@enumx NonlinearSafeTerminationReturnCode begin
    """
        NonlinearSafeTerminationReturnCode.Success

    Termination Condition was satisfied!
    """
    Success
    """
        NonlinearSafeTerminationReturnCode.Default

    Default Return Code. Used for type stability and conveys no additional information!
    """
    Default
    """
        NonlinearSafeTerminationReturnCode.PatienceTermination

    Terminate if there has been no improvement for the last `patience_steps`.
    """
    PatienceTermination
    """
        NonlinearSafeTerminationReturnCode.ProtectiveTermination

    Terminate if the objective value increased by this factor wrt initial objective or the
    value diverged.
    """
    ProtectiveTermination
    """
        NonlinearSafeTerminationReturnCode.Failure

    Termination Condition was not satisfied!
    """
    Failure
end

abstract type AbstractNonlinearTerminationMode end
abstract type AbstractSafeNonlinearTerminationMode <: AbstractNonlinearTerminationMode end
abstract type AbstractSafeBestNonlinearTerminationMode <:
              AbstractSafeNonlinearTerminationMode end

# TODO: Add a mode where the user can pass in custom termination criteria function

"""
    SteadyStateDiffEqTerminationMode <: AbstractNonlinearTerminationMode

Check if all values of the derivative is close to zero wrt both relative and absolute
tolerance.

The default used in SteadyStateDiffEq.jl! Not recommended for large problems, since the
convergence criteria is very strict and never reliably satisfied for most problems.
"""
struct SteadyStateDiffEqTerminationMode <: AbstractNonlinearTerminationMode end

"""
    SimpleNonlinearSolveTerminationMode <: AbstractNonlinearTerminationMode

Check if all values of the derivative is close to zero wrt both relative and absolute
tolerance. Or check that the value of the current and previous state is within the specified
tolerances.

The default used in SimpleNonlinearSolve.jl! Not recommended for large problems, since the
convergence criteria is very strict and never reliably satisfied for most problems.
"""
struct SimpleNonlinearSolveTerminationMode <: AbstractNonlinearTerminationMode end

@doc doc"""
    NormTerminationMode <: AbstractNonlinearTerminationMode

Terminates if
``\| \frac{\partial u}{\partial t} \| \leq reltol \times \| \frac{\partial u}{\partial t} + u \|``
or ``\| \frac{\partial u}{\partial t} \| \leq abstol``
"""
struct NormTerminationMode <: AbstractNonlinearTerminationMode end

@doc doc"""
    RelTerminationMode <: AbstractNonlinearTerminationMode

Terminates if
``all \left(| \frac{\partial u}{\partial t} | \leq reltol \times | u | \right)``.
"""
struct RelTerminationMode <: AbstractNonlinearTerminationMode end

@doc doc"""
    RelNormTerminationMode <: AbstractNonlinearTerminationMode

Terminates if
``\| \frac{\partial u}{\partial t} \| \leq reltol \times \| \frac{\partial u}{\partial t} + u \|``
"""
struct RelNormTerminationMode <: AbstractNonlinearTerminationMode end

@doc doc"""
    AbsTerminationMode <: AbstractNonlinearTerminationMode

Terminates if ``all \left( | \frac{\partial u}{\partial t} | \leq abstol \right)``.
"""
struct AbsTerminationMode <: AbstractNonlinearTerminationMode end

@doc doc"""
    AbsNormTerminationMode <: AbstractNonlinearTerminationMode

Terminates if ``\| \frac{\partial u}{\partial t} \| \leq abstol``.
"""
struct AbsNormTerminationMode <: AbstractNonlinearTerminationMode end

@doc doc"""
    RelSafeTerminationMode <: AbstractSafeNonlinearTerminationMode

Essentially [`RelNormTerminationMode`](@ref) + terminate if there has been no improvement
for the last `patience_steps` + terminate if the solution blows up (diverges).

## Constructor

```julia
RelSafeTerminationMode(; protective_threshold = nothing, patience_steps = 100,
                         patience_objective_multiplier = 3, min_max_factor = 1.3)
```
"""
Base.@kwdef struct RelSafeTerminationMode{T1, T2, T3} <:
                   AbstractSafeNonlinearTerminationMode
    protective_threshold::T1 = nothing
    patience_steps::Int = 100
    patience_objective_multiplier::T2 = 3
    min_max_factor::T3 = 1.3
end

@doc doc"""
    AbsSafeTerminationMode <: AbstractSafeNonlinearTerminationMode

Essentially [`AbsNormTerminationMode`](@ref) + terminate if there has been no improvement
for the last `patience_steps` + terminate if the solution blows up (diverges).

## Constructor

```julia
AbsSafeTerminationMode(; protective_threshold = nothing, patience_steps = 100,
                        patience_objective_multiplier = 3, min_max_factor = 1.3)
```
"""
Base.@kwdef struct AbsSafeTerminationMode{T1, T2, T3} <:
                   AbstractSafeNonlinearTerminationMode
    protective_threshold::T1 = nothing
    patience_steps::Int = 100
    patience_objective_multiplier::T2 = 3
    min_max_factor::T3 = 1.3
end

@doc doc"""
    RelSafeBestTerminationMode <: AbstractSafeBestNonlinearTerminationMode

Essentially [`RelSafeTerminationMode`](@ref), but caches the best solution found so far.

## Constructor

```julia
RelSafeBestTerminationMode(; protective_threshold = nothing, patience_steps = 100,
                            patience_objective_multiplier = 3, min_max_factor = 1.3)
```
"""
Base.@kwdef struct RelSafeBestTerminationMode{T1, T2, T3} <:
                   AbstractSafeBestNonlinearTerminationMode
    protective_threshold::T1 = nothing
    patience_steps::Int = 100
    patience_objective_multiplier::T2 = 3
    min_max_factor::T3 = 1.3
end

@doc doc"""
    AbsSafeBestTerminationMode <: AbstractSafeBestNonlinearTerminationMode

Essentially [`AbsSafeTerminationMode`](@ref), but caches the best solution found so far.

## Constructor

```julia
AbsSafeBestTerminationMode(; protective_threshold = nothing, patience_steps = 100,
                            patience_objective_multiplier = 3, min_max_factor = 1.3)
```
"""
Base.@kwdef struct AbsSafeBestTerminationMode{T1, T2, T3} <:
                   AbstractSafeBestNonlinearTerminationMode
    protective_threshold::T1 = nothing
    patience_steps::Int = 100
    patience_objective_multiplier::T2 = 3
    min_max_factor::T3 = 1.3
end

mutable struct NonlinearTerminationModeCache{uType, T,
    M <: AbstractNonlinearTerminationMode, I, OT, SV}
    u::uType
    retcode::NonlinearSafeTerminationReturnCode.T
    abstol::T
    reltol::T
    best_objective_value::T
    mode::M
    initial_objective::I
    objectives_trace::OT
    nsteps::Int
    saved_values::SV
end

get_termination_mode(cache::NonlinearTerminationModeCache) = cache.mode
get_abstol(cache::NonlinearTerminationModeCache) = cache.abstol
get_reltol(cache::NonlinearTerminationModeCache) = cache.reltol
get_saved_values(cache::NonlinearTerminationModeCache) = cache.saved_values

function __update_u!!(cache::NonlinearTerminationModeCache, u)
    cache.u === nothing && return
    if cache.u isa AbstractArray && ArrayInterface.can_setindex(cache.u)
        copyto!(cache.u, u)
    else
        cache.u = u
    end
end

__cvt_real(::Type{T}, ::Nothing) where {T} = nothing
__cvt_real(::Type{T}, x) where {T} = real(T(x))

_get_tolerance(η, ::Type{T}) where {T} = __cvt_real(T, η)
function _get_tolerance(::Nothing, ::Type{T}) where {T}
    η = real(oneunit(T)) * (eps(real(one(T))))^(4 // 5)
    return _get_tolerance(η, T)
end

function SciMLBase.init(du::Union{AbstractArray{T}, T}, u::Union{AbstractArray{T}, T},
        mode::AbstractNonlinearTerminationMode, saved_value_prototype...;
        abstol = nothing, reltol = nothing, kwargs...) where {T <: Number}
    abstol = _get_tolerance(abstol, T)
    reltol = _get_tolerance(reltol, T)
    TT = typeof(abstol)
    u_ = mode isa AbstractSafeBestNonlinearTerminationMode ?
         (ArrayInterface.can_setindex(u) ? copy(u) : u) : nothing
    if mode isa AbstractSafeNonlinearTerminationMode
        if mode isa AbsSafeTerminationMode || mode isa AbsSafeBestTerminationMode
            initial_objective = NONLINEARSOLVE_DEFAULT_NORM(du)
        else
            initial_objective = NONLINEARSOLVE_DEFAULT_NORM(du) /
                                (NONLINEARSOLVE_DEFAULT_NORM(du .+ u) + eps(TT))
        end
        objectives_trace = Vector{TT}(undef, mode.patience_steps)
        best_value = initial_objective
    else
        initial_objective = nothing
        objectives_trace = nothing
        best_value = __cvt_real(T, Inf)
    end

    length(saved_value_prototype) == 0 && (saved_value_prototype = nothing)

    return NonlinearTerminationModeCache{typeof(u_), TT, typeof(mode),
        typeof(initial_objective), typeof(objectives_trace),
        typeof(saved_value_prototype)}(u_, NonlinearSafeTerminationReturnCode.Default,
        abstol, reltol, best_value, mode, initial_objective, objectives_trace, 0,
        saved_value_prototype)
end

# This dispatch is needed based on how Terminating Callback works!
# This intentially drops the `abstol` and `reltol` arguments
function (cache::NonlinearTerminationModeCache)(integrator::AbstractODEIntegrator,
        abstol::Number, reltol::Number, min_t)
    retval = cache(cache.mode, get_du(integrator), integrator.u, integrator.uprev)
    (min_t === nothing || integrator.t ≥ min_t) && return retval
    return false
end
function (cache::NonlinearTerminationModeCache)(du, u, uprev, args...)
    return cache(cache.mode, du, u, uprev, args...)
end

function (cache::NonlinearTerminationModeCache)(mode::AbstractNonlinearTerminationMode, du,
        u, uprev, args...)
    return check_convergence(mode, du, u, uprev, cache.abstol, cache.reltol)
end

function (cache::NonlinearTerminationModeCache)(mode::AbstractSafeNonlinearTerminationMode,
        du, u, uprev, args...)
    if mode isa AbsSafeTerminationMode || mode isa AbsSafeBestTerminationMode
        objective = NONLINEARSOLVE_DEFAULT_NORM(du)
        criteria = cache.abstol
    else
        objective = NONLINEARSOLVE_DEFAULT_NORM(du) /
                    (NONLINEARSOLVE_DEFAULT_NORM(du .+ u) + eps(cache.abstol))
        criteria = cache.reltol
    end

    # Protective Break
    if isinf(objective) || isnan(objective)
        cache.retcode = NonlinearSafeTerminationReturnCode.ProtectiveTermination
        return true
    end
    ## By default we turn this off since it has the potential for false positives
    if cache.mode.protective_threshold !== nothing &&
       (objective > cache.initial_objective * cache.mode.protective_threshold * length(du))
        cache.retcode = NonlinearSafeTerminationReturnCode.ProtectiveTermination
        return true
    end

    # Check if best solution
    if mode isa AbstractSafeBestNonlinearTerminationMode &&
       objective < cache.best_objective_value
        cache.best_objective_value = objective
        __update_u!!(cache, u)
        if cache.saved_values !== nothing && length(args) ≥ 1
            cache.saved_values = args
        end
    end

    # Main Termination Condition
    if objective ≤ criteria
        cache.retcode = NonlinearSafeTerminationReturnCode.Success
        return true
    end

    # Terminate if there has been no improvement for the last `patience_steps`
    cache.nsteps += 1
    cache.nsteps == 1 && (cache.initial_objective = objective)
    cache.objectives_trace[mod1(cache.nsteps, length(cache.objectives_trace))] = objective

    if objective ≤ cache.mode.patience_objective_multiplier * criteria
        if cache.nsteps ≥ cache.mode.patience_steps
            if cache.nsteps < length(cache.objectives_trace)
                min_obj, max_obj = extrema(@view(cache.objectives_trace[1:(cache.nsteps)]))
            else
                min_obj, max_obj = extrema(cache.objectives_trace)
            end
            if min_obj < cache.mode.min_max_factor * max_obj
                cache.retcode = NonlinearSafeTerminationReturnCode.PatienceTermination
                return true
            end
        end
    end

    cache.retcode = NonlinearSafeTerminationReturnCode.Failure
    return false
end

const ZIPPABLE_TYPES = Union{Array, StaticArraysCore.StaticArray}

# Nonallocating version of `isapprox` if possible
function __nonlinearsolve_is_approx(x::ZIPPABLE_TYPES, y::ZIPPABLE_TYPES, abstol, reltol)
    length(x) != length(y) && return false
    # zip doesn't check lengths
    d = NONLINEARSOLVE_DEFAULT_NORM(((xᵢ, yᵢ),) -> xᵢ - yᵢ, zip(x, y))
    return d ≤ max(abstol,
        reltol * max(NONLINEARSOLVE_DEFAULT_NORM(x),
            NONLINEARSOLVE_DEFAULT_NORM(y)))
end
function __nonlinearsolve_is_approx(x, y, abstol, reltol)
    return isapprox(x, y; atol = abstol, rtol = reltol, norm = NONLINEARSOLVE_DEFAULT_NORM)
end

function check_convergence(::SteadyStateDiffEqTerminationMode, duₙ::ZIPPABLE_TYPES,
        uₙ::ZIPPABLE_TYPES, uₙ₋₁::ZIPPABLE_TYPES, abstol, reltol)
    return all(((x, y),) -> (abs(x) ≤ abstol) | (abs(x) ≤ reltol * abs(y)), zip(duₙ, uₙ))
end
function check_convergence(::SteadyStateDiffEqTerminationMode, duₙ, uₙ, uₙ₋₁, abstol,
        reltol)
    return all(@. (abs(duₙ) ≤ abstol) | (abs(duₙ) ≤ reltol * abs(uₙ)))
end

function check_convergence(::SimpleNonlinearSolveTerminationMode, duₙ::ZIPPABLE_TYPES,
        uₙ::ZIPPABLE_TYPES, uₙ₋₁::ZIPPABLE_TYPES, abstol, reltol)
    return all(((x, y),) -> (abs(x) ≤ abstol) | (abs(x) ≤ reltol * abs(y)), zip(duₙ, uₙ)) ||
           __nonlinearsolve_is_approx(uₙ, uₙ₋₁, abstol, reltol)  # isapprox allocates
end
function check_convergence(::SimpleNonlinearSolveTerminationMode, duₙ, uₙ, uₙ₋₁, abstol,
        reltol)
    return all(@. (abs(duₙ) ≤ abstol) | (abs(duₙ) ≤ reltol * abs(uₙ))) ||
           __nonlinearsolve_is_approx(uₙ, uₙ₋₁, abstol, reltol)  # isapprox allocates
end

function check_convergence(::NormTerminationMode, duₙ::ZIPPABLE_TYPES, uₙ::ZIPPABLE_TYPES,
        uₙ₋₁::ZIPPABLE_TYPES, abstol, reltol)
    du_norm = NONLINEARSOLVE_DEFAULT_NORM(duₙ)
    return du_norm ≤ abstol ||
           du_norm ≤ reltol * NONLINEARSOLVE_DEFAULT_NORM(((x, y),) -> x + y, zip(duₙ, uₙ))
end
function check_convergence(::NormTerminationMode, duₙ, uₙ, uₙ₋₁, abstol, reltol)
    du_norm = NONLINEARSOLVE_DEFAULT_NORM(duₙ)
    return du_norm ≤ abstol || du_norm ≤ reltol * NONLINEARSOLVE_DEFAULT_NORM(duₙ .+ uₙ)
end

function check_convergence(::RelTerminationMode, duₙ::ZIPPABLE_TYPES, uₙ::ZIPPABLE_TYPES,
        uₙ₋₁::ZIPPABLE_TYPES, abstol, reltol)
    return all(((x, y),) -> abs(x) ≤ reltol * abs(y), zip(duₙ, uₙ))
end
function check_convergence(::RelTerminationMode, duₙ, uₙ, uₙ₋₁, abstol, reltol)
    return all(@. abs(duₙ) ≤ reltol * abs(uₙ))
end

function check_convergence(::Union{RelNormTerminationMode, RelSafeTerminationMode,
            RelSafeBestTerminationMode}, duₙ::ZIPPABLE_TYPES, uₙ::ZIPPABLE_TYPES,
        uₙ₋₁::ZIPPABLE_TYPES, abstol, reltol)
    return NONLINEARSOLVE_DEFAULT_NORM(duₙ) ≤
           reltol * NONLINEARSOLVE_DEFAULT_NORM(((x, y),) -> x + y, zip(duₙ, uₙ))
end
function check_convergence(::Union{RelNormTerminationMode, RelSafeTerminationMode,
            RelSafeBestTerminationMode}, duₙ, uₙ, uₙ₋₁, abstol, reltol)
    return NONLINEARSOLVE_DEFAULT_NORM(duₙ) ≤
           reltol * NONLINEARSOLVE_DEFAULT_NORM(duₙ .+ uₙ)
end

function check_convergence(::AbsTerminationMode, duₙ, uₙ, uₙ₋₁, abstol, reltol)
    return all(x -> abs(x) ≤ abstol, duₙ)
end
function check_convergence(::Union{AbsNormTerminationMode, AbsSafeTerminationMode,
            AbsSafeBestTerminationMode}, duₙ, uₙ, uₙ₋₁, abstol, reltol)
    return NONLINEARSOLVE_DEFAULT_NORM(duₙ) ≤ abstol
end

# NOTE: Deprecate the following API eventually. This API leads to quite a bit of type
#       instability
@enumx NLSolveSafeTerminationReturnCode begin
    Success
    PatienceTermination
    ProtectiveTermination
    Failure
end

# SteadyStateDefault and NLSolveDefault are needed to be compatible with the existing
# termination conditions in NonlinearSolve and SteadyStateDiffEq
@enumx NLSolveTerminationMode begin
    SteadyStateDefault
    NLSolveDefault
    Norm
    Rel
    RelNorm
    Abs
    AbsNorm
    RelSafe
    RelSafeBest
    AbsSafe
    AbsSafeBest
end

struct NLSolveSafeTerminationOptions{T1, T2, T3}
    protective_threshold::T1
    patience_steps::Int
    patience_objective_multiplier::T2
    min_max_factor::T3
end

TruncatedStacktraces.@truncate_stacktrace NLSolveSafeTerminationOptions

mutable struct NLSolveSafeTerminationResult{T, uType}
    u::uType
    best_objective_value::T
    best_objective_value_iteration::Int
    return_code::NLSolveSafeTerminationReturnCode.T
end

function NLSolveSafeTerminationResult(u = nothing; best_objective_value = Inf64,
        best_objective_value_iteration = 0,
        return_code = NLSolveSafeTerminationReturnCode.Failure)
    u = u !== nothing ? copy(u) : u
    Base.depwarn("NLSolveSafeTerminationResult has been deprecated in favor of the new dispatch based termination conditions. Please use the new termination conditions API!",
        :NLSolveSafeTerminationResult)
    return NLSolveSafeTerminationResult{typeof(best_objective_value), typeof(u)}(u,
        best_objective_value, best_objective_value_iteration, return_code)
end

const BASIC_TERMINATION_MODES = (NLSolveTerminationMode.SteadyStateDefault,
    NLSolveTerminationMode.NLSolveDefault,
    NLSolveTerminationMode.Norm, NLSolveTerminationMode.Rel,
    NLSolveTerminationMode.RelNorm,
    NLSolveTerminationMode.Abs, NLSolveTerminationMode.AbsNorm)

const SAFE_TERMINATION_MODES = (NLSolveTerminationMode.RelSafe,
    NLSolveTerminationMode.RelSafeBest,
    NLSolveTerminationMode.AbsSafe,
    NLSolveTerminationMode.AbsSafeBest)

const SAFE_BEST_TERMINATION_MODES = (NLSolveTerminationMode.RelSafeBest,
    NLSolveTerminationMode.AbsSafeBest)

@doc doc"""
    NLSolveTerminationCondition(mode; abstol::T = 1e-8, reltol::T = 1e-6,
                                protective_threshold = 1e3, patience_steps::Int = 30,
                                patience_objective_multiplier = 3, min_max_factor = 1.3)

Define the termination criteria for the NonlinearProblem or SteadyStateProblem.

## Termination Conditions

#### Termination on Absolute Tolerance

  * `NLSolveTerminationMode.Abs`: Terminates if ``all \left( | \frac{\partial u}{\partial t} | \leq abstol \right)``
  * `NLSolveTerminationMode.AbsNorm`: Terminates if ``\| \frac{\partial u}{\partial t} \| \leq abstol``
  * `NLSolveTerminationMode.AbsSafe`: Essentially `abs_norm` + terminate if there has been no improvement for the last 30 steps + terminate if the solution blows up (diverges)
  * `NLSolveTerminationMode.AbsSafeBest`: Same as `NLSolveTerminationMode.AbsSafe` but uses the best solution found so far, i.e. deviates only if the solution has not converged

#### Termination on Relative Tolerance

  * `NLSolveTerminationMode.Rel`: Terminates if ``all \left(| \frac{\partial u}{\partial t} | \leq reltol \times | u | \right)``
  * `NLSolveTerminationMode.RelNorm`: Terminates if ``\| \frac{\partial u}{\partial t} \| \leq reltol \times \| \frac{\partial u}{\partial t} + u \|``
  * `NLSolveTerminationMode.RelSafe`: Essentially `rel_norm` + terminate if there has been no improvement for the last 30 steps + terminate if the solution blows up (diverges)
  * `NLSolveTerminationMode.RelSafeBest`: Same as `NLSolveTerminationMode.RelSafe` but uses the best solution found so far, i.e. deviates only if the solution has not converged

#### Termination using both Absolute and Relative Tolerances

  * `NLSolveTerminationMode.Norm`: Terminates if ``\| \frac{\partial u}{\partial t} \| \leq reltol \times \| \frac{\partial u}{\partial t} + u \|`` or ``\| \frac{\partial u}{\partial t} \| \leq abstol``
  * `NLSolveTerminationMode.SteadyStateDefault`: Check if all values of the derivative is close to zero wrt both relative and absolute tolerance. This is usable for small problems but doesn't scale well for neural networks.
  * `NLSolveTerminationMode.NLSolveDefault`: Check if all values of the derivative is close to zero wrt both relative and absolute tolerance. Or check that the value of the current and previous state is within the specified tolerances. This is usable for small problems but doesn't scale well for neural networks.

## General Arguments

  * `abstol`: Absolute Tolerance
  * `reltol`: Relative Tolerance

## Arguments specific to `*Safe*` modes

  * `protective_threshold`: If the objective value increased by this factor wrt initial objective terminate immediately.
  * `patience_steps`: If objective is within `patience_objective_multiplier` factor of the criteria and no improvement within `min_max_factor` has happened then terminate.

!!! warning
    This has been deprecated and will be removed in the next major release. Please use the new dispatch based termination conditions API.
"""
struct NLSolveTerminationCondition{mode, T,
    S <: Union{<:NLSolveSafeTerminationOptions, Nothing}}
    abstol::T
    reltol::T
    safe_termination_options::S
end

TruncatedStacktraces.@truncate_stacktrace NLSolveTerminationCondition 1

function Base.show(io::IO, s::NLSolveTerminationCondition{mode}) where {mode}
    print(io,
        "NLSolveTerminationCondition(mode = $(mode), abstol = $(s.abstol), reltol = $(s.reltol)")
    if mode ∈ SAFE_TERMINATION_MODES
        print(io, ", safe_termination_options = ", s.safe_termination_options, ")")
    else
        print(io, ")")
    end
end

get_termination_mode(::NLSolveTerminationCondition{mode}) where {mode} = mode

# Don't specify `mode` since the defaults would depend on the package
function NLSolveTerminationCondition(mode; abstol::T = 1e-8, reltol::T = 1e-6,
        protective_threshold = 1e3, patience_steps::Int = 30,
        patience_objective_multiplier = 3,
        min_max_factor = 1.3) where {T}
    Base.depwarn("NLSolveTerminationCondition has been deprecated in favor of the new dispatch based termination conditions. Please use the new termination conditions API!",
        :NLSolveTerminationCondition)
    @assert mode ∈ instances(NLSolveTerminationMode.T)
    options = if mode ∈ SAFE_TERMINATION_MODES
        NLSolveSafeTerminationOptions(protective_threshold, patience_steps,
            patience_objective_multiplier, min_max_factor)
    else
        nothing
    end
    return NLSolveTerminationCondition{mode, T, typeof(options)}(abstol, reltol, options)
end

function (cond::NLSolveTerminationCondition)(storage::Union{
        NLSolveSafeTerminationResult,
        Nothing,
    })
    mode = get_termination_mode(cond)
    # We need both the dispatches to support solvers that don't use the integrator
    # interface like SimpleNonlinearSolve
    if mode in BASIC_TERMINATION_MODES
        function _termination_condition_closure_basic(integrator, abstol, reltol, min_t)
            return _termination_condition_closure_basic(get_du(integrator), integrator.u,
                integrator.uprev, abstol, reltol)
        end
        function _termination_condition_closure_basic(du, u, uprev, abstol, reltol)
            return _has_converged(du, u, uprev, cond, abstol, reltol)
        end
        return _termination_condition_closure_basic
    else
        mode ∈ SAFE_BEST_TERMINATION_MODES && @assert storage !== nothing
        nstep::Int = 0

        function _termination_condition_closure_safe(integrator, abstol, reltol, min_t)
            return _termination_condition_closure_safe(get_du(integrator), integrator.u,
                integrator.uprev, abstol, reltol)
        end
        @inbounds function _termination_condition_closure_safe(du, u, uprev, abstol, reltol)
            aType = typeof(abstol)
            protective_threshold = aType(cond.safe_termination_options.protective_threshold)
            objective_values = aType[]
            patience_objective_multiplier = cond.safe_termination_options.patience_objective_multiplier

            if mode ∈ SAFE_BEST_TERMINATION_MODES
                storage.best_objective_value = aType(Inf)
                storage.best_objective_value_iteration = 0
            end

            if mode ∈ SAFE_BEST_TERMINATION_MODES
                objective = NONLINEARSOLVE_DEFAULT_NORM(du)
                criteria = abstol
            else
                objective = NONLINEARSOLVE_DEFAULT_NORM(du) /
                            (NONLINEARSOLVE_DEFAULT_NORM(du .+ u) + eps(aType))
                criteria = reltol
            end

            if mode ∈ SAFE_BEST_TERMINATION_MODES
                if objective < storage.best_objective_value
                    storage.best_objective_value = objective
                    storage.best_objective_value_iteration = nstep + 1
                    if storage.u !== nothing
                        storage.u .= u
                    end
                end
            end

            # Main Termination Criteria
            if objective ≤ criteria
                storage.return_code = NLSolveSafeTerminationReturnCode.Success
                return true
            end

            # Terminate if there has been no improvement for the last `patience_steps`
            nstep += 1
            push!(objective_values, objective)

            if objective ≤ typeof(criteria)(patience_objective_multiplier) * criteria
                if nstep ≥ cond.safe_termination_options.patience_steps
                    last_k_values = objective_values[max(1,
                        length(objective_values) -
                        cond.safe_termination_options.patience_steps):end]
                    if maximum(last_k_values) <
                       typeof(criteria)(cond.safe_termination_options.min_max_factor) *
                       minimum(last_k_values)
                        storage.return_code = NLSolveSafeTerminationReturnCode.PatienceTermination
                        return true
                    end
                end
            end

            # Protective break
            if objective ≥ objective_values[1] * protective_threshold * length(du)
                storage.return_code = NLSolveSafeTerminationReturnCode.ProtectiveTermination
                return true
            end

            storage.return_code = NLSolveSafeTerminationReturnCode.Failure
            return false
        end
        return _termination_condition_closure_safe
    end
end

# Convergence Criteria
@inline function _has_converged(du, u, uprev, cond::NLSolveTerminationCondition{mode},
        abstol = cond.abstol, reltol = cond.reltol) where {mode}
    return _has_converged(du, u, uprev, mode, abstol, reltol)
end

@inline @inbounds function _has_converged(du, u, uprev, mode, abstol, reltol)
    if mode == NLSolveTerminationMode.Norm
        du_norm = NONLINEARSOLVE_DEFAULT_NORM(du)
        return du_norm ≤ abstol || du_norm ≤ reltol * NONLINEARSOLVE_DEFAULT_NORM(du + u)
    elseif mode == NLSolveTerminationMode.Rel
        return all(abs.(du) .≤ reltol .* abs.(u))
    elseif mode ∈ (NLSolveTerminationMode.RelNorm, NLSolveTerminationMode.RelSafe,
        NLSolveTerminationMode.RelSafeBest)
        return NONLINEARSOLVE_DEFAULT_NORM(du) ≤
               reltol * NONLINEARSOLVE_DEFAULT_NORM(du .+ u)
    elseif mode == NLSolveTerminationMode.Abs
        return all(abs.(du) .≤ abstol)
    elseif mode ∈ (NLSolveTerminationMode.AbsNorm, NLSolveTerminationMode.AbsSafe,
        NLSolveTerminationMode.AbsSafeBest)
        return NONLINEARSOLVE_DEFAULT_NORM(du) ≤ abstol
    elseif mode == NLSolveTerminationMode.SteadyStateDefault
        return all((abs.(du) .≤ abstol) .| (abs.(du) .≤ reltol .* abs.(u)))
    elseif mode == NLSolveTerminationMode.NLSolveDefault
        atol, rtol = abstol, reltol
        return all((abs.(du) .≤ abstol) .| (abs.(du) .≤ reltol .* abs.(u))) ||
               isapprox(u, uprev; atol, rtol)
    else
        throw(ArgumentError("Unknown termination mode: $mode"))
    end
end
