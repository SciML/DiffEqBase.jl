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

"""
struct NLSolveTerminationCondition{mode, T,
                                   S <: Union{<:NLSolveSafeTerminationOptions, Nothing}}
    abstol::T
    reltol::T
    safe_termination_options::S
end

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
    @assert mode ∈ instances(NLSolveTerminationMode.T)
    options = if mode ∈ SAFE_TERMINATION_MODES
        NLSolveSafeTerminationOptions(protective_threshold, patience_steps,
                                      patience_objective_multiplier, min_max_factor)
    else
        nothing
    end
    return NLSolveTerminationCondition{mode, T, typeof(options)}(abstol, reltol, options)
end

function (cond::NLSolveTerminationCondition)(storage::Union{<:AbstractDict, Nothing})
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
                storage[:best_objective_value] = aType(Inf)
                storage[:best_objective_value_iteration] = 0
            end

            if mode ∈ SAFE_BEST_TERMINATION_MODES
                objective = norm(du)
                criteria = abstol
            else
                objective = norm(du) / (norm(du .+ u) + eps(aType))
                criteria = reltol
            end

            if mode ∈ SAFE_BEST_TERMINATION_MODES
                if objective < storage[:best_objective_value]
                    storage[:best_objective_value] = objective
                    storage[:best_objective_value_iteration] = nstep + 1
                end
            end

            # Main Termination Criteria
            if objective <= criteria
                storage[:return_code] = NLSolveSafeTerminationReturnCode.Success
                return true
            end

            # Terminate if there has been no improvement for the last `patience_steps`
            nstep += 1
            push!(objective_values, objective)

            if objective <= typeof(criteria)(patience_objective_multiplier) * criteria
                if nstep >= cond.safe_termination_options.patience_steps
                    last_k_values = objective_values[max(1,
                                                         length(objective_values) -
                                                         cond.safe_termination_options.patience_steps):end]
                    if maximum(last_k_values) <
                       typeof(criteria)(cond.safe_termination_options.min_max_factor) *
                       minimum(last_k_values)
                        storage[:return_code] = NLSolveSafeTerminationReturnCode.PatienceTermination
                        return true
                    end
                end
            end

            # Protective break
            if objective >= objective_values[1] * protective_threshold * length(du)
                storage[:return_code] = NLSolveSafeTerminationReturnCode.ProtectiveTermination
                return true
            end

            storage[:return_code] = NLSolveSafeTerminationReturnCode.Failure
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
        du_norm = norm(du)
        return du_norm <= abstol || du_norm <= reltol * norm(du + u)
    elseif mode == NLSolveTerminationMode.Rel
        return all(abs.(du) .<= reltol .* abs.(u))
    elseif mode ∈ (NLSolveTerminationMode.RelNorm, NLSolveTerminationMode.RelSafe,
            NLSolveTerminationMode.RelSafeBest)
        return norm(du) <= reltol * norm(du .+ u)
    elseif mode == NLSolveTerminationMode.Abs
        return all(abs.(du) .<= abstol)
    elseif mode ∈ (NLSolveTerminationMode.AbsNorm, NLSolveTerminationMode.AbsSafe,
            NLSolveTerminationMode.AbsSafeBest)
        return norm(du) <= abstol
    elseif mode == NLSolveTerminationMode.SteadyStateDefault
        return all((abs.(du) .<= abstol) .| (abs.(du) .<= reltol .* abs.(u)))
    elseif mode == NLSolveTerminationMode.NLSolveDefault
        atol, rtol = abstol, reltol
        return all((abs.(du) .<= abstol) .| (abs.(du) .<= reltol .* abs.(u))) ||
               isapprox(u, uprev; atol, rtol)
    else
        throw(ArgumentError("Unknown termination mode: $mode"))
    end
end
