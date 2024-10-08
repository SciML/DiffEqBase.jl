module DiffEqBase
if isdefined(Base, :Experimental) &&
   isdefined(Base.Experimental, Symbol("@max_methods"))
    @eval Base.Experimental.@max_methods 1
end

import PrecompileTools

import FastPower
@deprecate fastpow(x,y) FastPower.fastpower(x,y)

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

using ForwardDiff

using EnumX

using Markdown

using ConcreteStructs: @concrete
using FastClosures: @closure

# Could be made optional/glue
import PreallocationTools

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
                 AbstractSDDEAlgorithm, AbstractRODEAlgorithm,
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
                 solution_new_tslocation, plot_indices,
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

import SciMLBase: AbstractDiffEqLinearOperator # deprecation path

import SciMLStructures

import Tricks

using Reexport
Reexport.@reexport using SciMLBase

SciMLBase.isfunctionwrapper(x::FunctionWrapper) = true

"""
$(TYPEDEF)
"""
abstract type Tableau end

"""
$(TYPEDEF)
"""
abstract type ODERKTableau <: Tableau end

"""
$(TYPEDEF)
"""
abstract type DECostFunction end

import SciMLBase: Void, unwrapped_f

include("utils.jl")
include("stats.jl")
include("calculate_residuals.jl")
include("tableaus.jl")
include("internal_falsi.jl")
include("internal_itp.jl")

include("callbacks.jl")
include("common_defaults.jl")
include("solve.jl")
include("internal_euler.jl")
include("forwarddiff.jl")
include("termination_conditions_deprecated.jl")  # TODO: remove in the next major release
include("termination_conditions.jl")
include("norecompile.jl")
include("integrator_accessors.jl")

# This is only used for oop stiff solvers
default_factorize(A) = lu(A; check = false)

if isdefined(SciMLBase, :AbstractParameterizedFunction)
    import SciMLBase: AbstractParameterizedFunction
else
    """
    $(TYPEDEF)
    """
    abstract type AbstractParameterizedFunction{iip} <: AbstractODEFunction{iip} end
end

"""
$(TYPEDEF)
"""
struct ConvergenceSetup{P, C}
    probs::P
    convergence_axis::C
end

export initialize!, finalize!

export SensitivityADPassThrough

export SteadyStateDiffEqTerminationMode, SimpleNonlinearSolveTerminationMode,
       NormTerminationMode, RelTerminationMode, RelNormTerminationMode, AbsTerminationMode,
       AbsNormTerminationMode, RelSafeTerminationMode, AbsSafeTerminationMode,
       RelSafeBestTerminationMode, AbsSafeBestTerminationMode
# Deprecated API
export NLSolveTerminationMode,
       NLSolveSafeTerminationOptions, NLSolveTerminationCondition,
       NLSolveSafeTerminationResult

export KeywordArgError, KeywordArgWarn, KeywordArgSilent

end # module
