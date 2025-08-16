module DiffEqBase
if isdefined(Base, :Experimental) &&
   isdefined(Base.Experimental, Symbol("@max_methods"))
    @eval Base.Experimental.@max_methods 1
end

import PrecompileTools

import FastPower
@deprecate fastpow(x, y) FastPower.fastpower(x, y)

using ArrayInterface

using StaticArraysCore # data arrays

using LinearAlgebra, Printf

using DocStringExtensions

using FunctionWrappers: FunctionWrapper

using MuladdMacro, Parameters

using Statistics

using FastBroadcast: @.., True, False

using Static: reduce_tup

import RecursiveArrayTools
import TruncatedStacktraces

using Setfield

using EnumX

using Markdown

using ConcreteStructs: @concrete
using FastClosures: @closure

import FunctionWrappersWrappers

using SciMLBase

using SciMLOperators: AbstractSciMLOperator, AbstractSciMLScalarOperator

using SciMLBase: @def, DEIntegrator, AbstractDEProblem,
                 AbstractDiffEqInterpolation,
                 DECallback, AbstractDEOptions, DECache, AbstractContinuousCallback,
                 AbstractDiscreteCallback, AbstractLinearProblem,
                 AbstractNonlinearProblem,
                 AbstractOptimizationProblem, AbstractSteadyStateProblem,
                 AbstractJumpProblem,
                 AbstractNoiseProblem, AbstractEnsembleProblem,
                 AbstractDynamicalODEProblem,
                 AbstractDEAlgorithm, StandardODEProblem, AbstractIntegralProblem,
                 AbstractSensitivityAlgorithm, AbstractODEAlgorithm,
                 AbstractSDEAlgorithm, AbstractDDEAlgorithm, AbstractDAEAlgorithm,
                 AbstractSDDEAlgorithm, AbstractRODEAlgorithm, AbstractBVPAlgorithm,
                 DAEInitializationAlgorithm,
                 AbstractSteadyStateAlgorithm, AbstractODEProblem,
                 AbstractDiscreteProblem, AbstractNonlinearAlgorithm,
                 AbstractSDEProblem, AbstractRODEProblem, AbstractDDEProblem,
                 AbstractDAEProblem, AbstractSDDEProblem, AbstractBVProblem,
                 AbstractTimeseriesSolution, AbstractNoTimeSolution, numargs,
                 AbstractODEFunction, AbstractSDEFunction, AbstractRODEFunction,
                 AbstractDDEFunction, AbstractSDDEFunction, AbstractDAEFunction,
                 AbstractNonlinearFunction, AbstractEnsembleSolution,
                 AbstractODESolution, AbstractRODESolution, AbstractDAESolution,
                 AbstractDDESolution,
                 EnsembleAlgorithm, EnsembleSolution, EnsembleSummary,
                 NonlinearSolution,
                 TimeGradientWrapper, TimeDerivativeWrapper, UDerivativeWrapper,
                 UJacobianWrapper, ParamJacobianWrapper, JacobianWrapper,
                 check_error!, has_jac, has_tgrad, has_Wfact, has_Wfact_t, has_paramjac,
                 AbstractODEIntegrator, AbstractSDEIntegrator, AbstractRODEIntegrator,
                 AbstractDDEIntegrator, AbstractSDDEIntegrator,
                 AbstractDAEIntegrator, unwrap_cache, has_reinit, reinit!,
                 postamble!, last_step_failed, islinear, has_stats,
                 initialize_dae!, build_solution, solution_new_retcode,
                 solution_new_tslocation, plot_indices, NonlinearAliasSpecifier,
                 NullParameters, isinplace, AbstractADType, AbstractDiscretization,
                 DISCRETE_OUTOFPLACE_DEFAULT, DISCRETE_INPLACE_DEFAULT,
                 has_analytic, calculate_solution_errors!, AbstractNoiseProcess,
                 has_colorvec, parameterless_type, undefined_exports,
                 is_diagonal_noise, AbstractDiffEqFunction, sensitivity_solution,
                 interp_summary, AbstractHistoryFunction, LinearInterpolation,
                 ConstantInterpolation, HermiteInterpolation, SensitivityInterpolation,
                 NoAD, @add_kwonly,
                 calculate_ensemble_errors, DEFAULT_UPDATE_FUNC, isconstant,
                 DEFAULT_REDUCTION, isautodifferentiable,
                 isadaptive, isdiscrete, has_syms, AbstractAnalyticalSolution,
                 RECOMPILE_BY_DEFAULT, wrap_sol, has_destats

import SciMLBase: solve, init, step!, solve!, __init, __solve, update_coefficients!,
                  update_coefficients, isadaptive, wrapfun_oop, wrapfun_iip,
                  unwrap_fw, promote_tspan, set_u!, set_t!, set_ut!

import SciMLStructures

using Reexport
Reexport.@reexport using SciMLBase

SciMLBase.isfunctionwrapper(x::FunctionWrapper) = true

import SymbolicIndexingInterface as SII

## Extension Functions

eltypedual(x) = false
promote_u0(::Nothing, p, t0) = nothing
isdualtype(::Type{T}) where {T} = false

## Types

"""
    Tableau

Abstract type for Butcher tableaus used in Runge-Kutta methods.

Tableaus define the coefficients for multi-stage integration methods,
including the `a` (coupling), `b` (weights), and `c` (nodes) coefficients.

# Subtypes
- `ODERKTableau`: Tableaus specifically for explicit Runge-Kutta ODE methods

# See Also
- [`ODERKTableau`](@ref)
"""
abstract type Tableau end

"""
    ODERKTableau <: Tableau

Abstract type for Butcher tableaus specific to explicit Runge-Kutta ODE methods.

These tableaus are used to define the coefficients for explicit RK methods,
where the coupling matrix `a` is strictly lower triangular.

# Fields (typically implemented by concrete types)
- `a`: Coupling coefficients matrix (lower triangular)
- `b`: Weight coefficients for final step
- `c`: Time nodes for intermediate stages

# See Also
- [`Tableau`](@ref)
"""
abstract type ODERKTableau <: Tableau end

"""
    DECostFunction

Abstract type for cost/objective functions used in optimization-based
differential equation solving methods.

Cost functions define the objective to minimize when solving DEs using
optimization approaches, such as in collocation methods or parameter estimation.

# Implementation
Concrete subtypes should define methods for evaluating the cost and its gradients.
"""
abstract type DECostFunction end

import SciMLBase: Void, unwrapped_f

include("utils.jl")
include("stats.jl")
include("calculate_residuals.jl")
include("tableaus.jl")
include("internal_itp.jl")

include("callbacks.jl")
include("common_defaults.jl")
include("solve.jl")
include("internal_euler.jl")
include("norecompile.jl")
include("integrator_accessors.jl")

# This is only used for oop stiff solvers
default_factorize(A) = lu(A; check = false)

if isdefined(SciMLBase, :AbstractParameterizedFunction)
    import SciMLBase: AbstractParameterizedFunction
else
    """
        AbstractParameterizedFunction{iip} <: AbstractODEFunction{iip}

    Abstract type for parameterized ODE functions that depend on additional parameters.

    Parameterized functions allow for ODEs where the dynamics explicitly depend on
    parameters that may be varied, fitted, or optimized.

    # Type Parameters
    - `iip`: Whether the function is in-place (true) or out-of-place (false)

    # See Also
    - [`AbstractODEFunction`](@ref)
    """
    abstract type AbstractParameterizedFunction{iip} <: AbstractODEFunction{iip} end
end

"""
    ConvergenceSetup{P, C}

Configuration for convergence analysis of differential equation solvers.

Used to test and verify the convergence properties of numerical methods by
solving a set of problems with varying discretization parameters.

# Fields
- `probs::P`: Collection of test problems with different discretizations
- `convergence_axis::C`: The axis along which convergence is measured (e.g., time step sizes)

# Usage
Typically used internally by convergence testing utilities to verify that
a solver achieves its expected order of accuracy.

# Example
```julia
probs = [ODEProblem(f, u0, tspan; dt=dt) for dt in [0.1, 0.05, 0.025, 0.0125]]
setup = ConvergenceSetup(probs, [0.1, 0.05, 0.025, 0.0125])
```
"""
struct ConvergenceSetup{P, C}
    probs::P
    convergence_axis::C
end

export initialize!, finalize!

export SensitivityADPassThrough

export KeywordArgError, KeywordArgWarn, KeywordArgSilent

end # module
