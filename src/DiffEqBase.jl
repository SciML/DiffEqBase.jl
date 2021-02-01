module DiffEqBase

using Requires, IterativeSolvers, RecursiveFactorization, ArrayInterface

using StaticArrays # data arrays

using LinearAlgebra, Printf

using DocStringExtensions

using FunctionWrappers: FunctionWrapper

using MuladdMacro, Parameters

using NonlinearSolve

using Reexport

using Statistics

import ZygoteRules, ChainRulesCore
import LabelledArrays
import RecursiveArrayTools

@reexport using SciMLBase

using SciMLBase: @def, DEIntegrator, DEProblem, AbstractDiffEqOperator,
                 AbstractDiffEqLinearOperator, AbstractDiffEqInterpolation,
                 DECallback, AbstractDEOptions, DECache, AbstractContinuousCallback,
                 AbstractDiscreteCallback, AbstractLinearProblem, AbstractNonlinearProblem,
                 AbstractOptimizationProblem, AbstractQuadratureProblem,
                 AbstractSteadyStateProblem, AbstractJumpProblem,
                 AbstractNoiseProblem, AbstractEnsembleProblem, AbstractDynamicalODEProblem,
                 DEAlgorithm,
                 AbstractSensitivityAlgorithm, AbstractODEAlgorithm,
                 AbstractSDEAlgorithm, AbstractDDEAlgorithm, AbstractDAEAlgorithm,
                 AbstractSDDEAlgorithm, AbstractRODEAlgorithm, DAEInitializationAlgorithm,
                 AbstractSteadyStateAlgorithm, AbstractQuadratureAlgorithm,
                 AbstractODEProblem, AbstractDiscreteProblem,
                 AbstractSDEProblem, AbstractRODEProblem, AbstractDDEProblem,
                 AbstractDAEProblem, AbstractSDDEProblem, AbstractBVProblem,
                 AbstractTimeseriesSolution, AbstractNoTimeSolution, numargs,
                 AbstractODEFunction, AbstractSDEFunction, AbstractRODEFunction,
                 AbstractDDEFunction, AbstractSDDEFunction, AbstractDAEFunction,
                 AbstractNonlinearFunction, AbstractEnsembleSolution,
                 AbstractODESolution, AbstractRODESolution, AbstractDAESolution,
                 EnsembleAlgorithm, EnsembleSolution, EnsembleSummary,
                 QuadratureSolution,
                 TimeGradientWrapper, TimeDerivativeWrapper, UDerivativeWrapper,
                 UJacobianWrapper, ParamJacobianWrapper, JacobianWrapper,
                 check_error!, has_jac, has_tgrad, has_Wfact, has_Wfact_t, has_paramjac,
                 AbstractODEIntegrator, AbstractSDEIntegrator, AbstractRODEIntegrator,
                 AbstractDDEIntegrator, AbstractSDDEIntegrator,
                 AbstractDAEIntegrator, unwrap_cache, has_reinit, reinit!,
                 postamble!, last_step_failed, islinear, has_destats,
                 initialize_dae!, build_solution, solution_new_retcode,
                 solution_new_tslocation, plot_indices,
                 NullParameters, isinplace, AbstractADType, AbstractDiscretization,
                 DISCRETE_OUTOFPLACE_DEFAULT, DISCRETE_INPLACE_DEFAULT,
                 has_analytic, calculate_solution_errors!, AbstractNoiseProcess,
                 has_colorvec, parameterless_type, undefined_exports,
                 is_diagonal_noise, AbstractDiffEqFunction, sensitivity_solution,
                 interp_summary, AbstractHistoryFunction, LinearInterpolation,
                 ConstantInterpolation, HermiteInterpolation, NoAD, @add_kwonly,
                 calculate_ensemble_errors, DEFAULT_UPDATE_FUNC, isconstant,
                 DEFAULT_REDUCTION, isautodifferentiable,
                 isadaptive, isdiscrete

import SciMLBase: solve, init, solve!, __init, __solve, update_coefficients!,
                  update_coefficients, isadaptive, wrapfun_oop, wrapfun_iip,
                  unwrap_fw, promote_tspan, set_u!, set_t!, set_ut!

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

"""
$(TYPEDEF)
"""
abstract type DEDataArray{T,N} <: AbstractArray{T,N} end
const DEDataVector{T} = DEDataArray{T,1}
const DEDataMatrix{T} = DEDataArray{T,2}

include("utils.jl")
include("fastpow.jl")
include("diffeqfastbc.jl")
include("destats.jl")
include("calculate_residuals.jl")
include("tableaus.jl")

include("nlsolve/type.jl")
include("nlsolve/newton.jl")
include("nlsolve/functional.jl")
include("nlsolve/utils.jl")
include("callbacks.jl")
include("linear_nonlinear.jl")
include("common_defaults.jl")
include("data_array.jl")
include("solve.jl")
include("internal_euler.jl")
include("init.jl")
include("zygote.jl")

"""
$(TYPEDEF)
"""
abstract type AbstractParameterizedFunction{iip} <: AbstractODEFunction{iip} end

"""
$(TYPEDEF)
"""
struct ConvergenceSetup{P,C}
    probs::P
    convergence_axis::C
end

const AbstractMonteCarloProblem = AbstractEnsembleProblem
const AbstractMonteCarloSolution = AbstractEnsembleSolution
const MonteCarloAlgorithm = EnsembleAlgorithm
const MonteCarloProblem = EnsembleProblem
const MonteCarloSolution = EnsembleSolution
const MonteCarloSummary = EnsembleSummary
@deprecate MonteCarloProblem(args...) EnsembleProblem(args...)
@deprecate MonteCarloSolution(args...) EnsembleSolution(args...)
@deprecate MonteCarloSummary(args...) EnsembleSummary(args...)
@deprecate calculate_monte_errors(args...;kwargs...) calculate_ensemble_errors(args...;kwargs...)

export DEDataArray, DEDataVector, DEDataMatrix

export ContinuousCallback, DiscreteCallback, CallbackSet, VectorContinuousCallback

export initialize!, finalize!

export LinSolveFactorize, LinSolveGPUFactorize, DefaultLinSolve, DEFAULT_LINSOLVE,
       LinSolveGMRES, LinSolveCG, LinSolveBiCGStabl, LinSolveChebyshev,
       LinSolveMINRES, LinSolveIterativeSolvers

export NLNewton, NLFunctional, NLAnderson

export SensitivityADPassThrough

end # module
