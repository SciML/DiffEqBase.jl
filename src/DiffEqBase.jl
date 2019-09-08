module DiffEqBase

using RecipesBase, RecursiveArrayTools, Compat,
      Requires, TableTraits, IteratorInterfaceExtensions, TreeViews,
      IterativeSolvers, RecursiveFactorization, Distributed, ArrayInterface

using Roots # callbacks

using StaticArrays # data arrays

using LinearAlgebra, Statistics, Printf

using DocStringExtensions

using FunctionWrappers: FunctionWrapper

using MuladdMacro, Parameters

# Problems
"""
$(TYPEDEF)

Base type for all DifferentialEquations.jl problems. Concrete subtypes of
`DEProblem` contain the necessary information to fully define a differential
equation of the corresponding type.
"""
abstract type DEProblem end

"""
$(TYPEDEF)
"""
abstract type DEElement end

"""
$(TYPEDEF)
"""
abstract type DESensitivity end

"""
$(TYPEDEF)

Base for types which define steady state problems for ODE systems.
"""
abstract type AbstractSteadyStateProblem{uType,isinplace} <: DEProblem end

"""
$(TYPEDEF)
"""
abstract type AbstractNoiseProblem <: DEProblem end

"""
$(TYPEDEF)

Base for types which define ODE problems.
"""
abstract type AbstractODEProblem{uType,tType,isinplace} <: DEProblem end

"""
$(TYPEDEF)

Base for types which define discrete problems.
"""
abstract type AbstractDiscreteProblem{uType,tType,isinplace} <:
                      AbstractODEProblem{uType,tType,isinplace} end

"""
$(TYPEDEF)
"""
abstract type AbstractAnalyticalProblem{uType,tType,isinplace} <:
                      AbstractODEProblem{uType,tType,isinplace} end

"""
$(TYPEDEF)

Base for types which define RODE problems.
"""
abstract type AbstractRODEProblem{uType,tType,isinplace,ND} <: DEProblem end

"""
$(TYPEDEF)

Base for types which define SDE problems.
"""
abstract type AbstractSDEProblem{uType,tType,isinplace,ND} <:
                      AbstractRODEProblem{uType,tType,isinplace,ND} end

"""
$(TYPEDEF)

Base for types which define DAE problems.
"""
abstract type AbstractDAEProblem{uType,duType,tType,isinplace} <: DEProblem end

"""
$(TYPEDEF)

Base for types which define DDE problems.
"""
abstract type AbstractDDEProblem{uType,tType,lType,isinplace} <: DEProblem end

"""
$(TYPEDEF)
"""
abstract type AbstractConstantLagDDEProblem{uType,tType,lType,isinplace} <:
                      AbstractDDEProblem{uType,tType,lType,isinplace} end

"""
$(TYPEDEF)
"""
abstract type AbstractSecondOrderODEProblem{uType,tType,isinplace} <: AbstractODEProblem{uType,tType,isinplace} end

"""
$(TYPEDEF)

Base for types which define BVP problems.
"""
abstract type AbstractBVProblem{uType,tType,isinplace} <: AbstractODEProblem{uType,tType,isinplace} end

"""
$(TYPEDEF)

Base for types which define jump problems.
"""
abstract type AbstractJumpProblem{P,J} <: DiffEqBase.DEProblem end

# Algorithms
"""
$(TYPEDEF)
"""
abstract type DEAlgorithm end

"""
$(TYPEDEF)
"""
abstract type AbstractSteadyStateAlgorithm <: DEAlgorithm end

"""
$(TYPEDEF)
"""
abstract type AbstractODEAlgorithm <: DEAlgorithm end

"""
$(TYPEDEF)
"""
abstract type AbstractSecondOrderODEAlgorithm <: DEAlgorithm end

"""
$(TYPEDEF)
"""
abstract type AbstractRODEAlgorithm <: DEAlgorithm end

"""
$(TYPEDEF)
"""
abstract type AbstractSDEAlgorithm <: DEAlgorithm end

"""
$(TYPEDEF)
"""
abstract type AbstractDAEAlgorithm <: DEAlgorithm end

"""
$(TYPEDEF)
"""
abstract type AbstractDDEAlgorithm <: DEAlgorithm end

"""
$(TYPEDEF)
"""
abstract type EnsembleAlgorithm <: DiffEqBase.DEAlgorithm end

# Monte Carlo Simulations
"""
$(TYPEDEF)
"""
abstract type AbstractEnsembleProblem <: DEProblem end

"""
$(TYPEDEF)
"""
abstract type AbstractEnsembleEstimator <: DEProblem end

export EnsembleProblem
export EnsembleSolution, EnsembleTestSolution, EnsembleSummary

"""
$(TYPEDEF)
"""
abstract type AbstractDiffEqInterpolation <: Function end

"""
$(TYPEDEF)
"""
abstract type AbstractDEOptions end

"""
$(TYPEDEF)
"""
abstract type DECache end

"""
$(TYPEDEF)
"""
abstract type DECallback end

"""
$(TYPEDEF)
"""
abstract type AbstractContinuousCallback <: DECallback end

"""
$(TYPEDEF)
"""
abstract type AbstractDiscreteCallback <: DECallback end

"""
$(TYPEDEF)
"""
abstract type DEDataArray{T,N} <: AbstractArray{T,N} end
const DEDataVector{T} = DEDataArray{T,1}
const DEDataMatrix{T} = DEDataArray{T,2}

# Integrators
"""
$(TYPEDEF)
"""
abstract type DEIntegrator{Alg, IIP, U, T} end

"""
$(TYPEDEF)
"""
abstract type AbstractSteadyStateIntegrator{Alg, IIP, U} <: DEIntegrator{Alg, IIP, U, Nothing} end

"""
$(TYPEDEF)
"""
abstract type AbstractODEIntegrator{Alg, IIP, U, T} <: DEIntegrator{Alg, IIP, U, T} end

"""
$(TYPEDEF)
"""
abstract type AbstractSecondOrderODEIntegrator{Alg, IIP, U, T} <: DEIntegrator{Alg, IIP, U, T} end

"""
$(TYPEDEF)
"""
abstract type AbstractRODEIntegrator{Alg, IIP, U, T} <: DEIntegrator{Alg, IIP, U, T} end

"""
$(TYPEDEF)
"""
abstract type AbstractSDEIntegrator{Alg, IIP, U, T} <: DEIntegrator{Alg, IIP, U, T} end

"""
$(TYPEDEF)
"""
abstract type AbstractDDEIntegrator{Alg, IIP, U, T} <: DEIntegrator{Alg, IIP, U, T} end

"""
$(TYPEDEF)
"""
abstract type AbstractDAEIntegrator{Alg, IIP, U, T} <: DEIntegrator{Alg, IIP, U, T} end

# Solutions
"""
$(TYPEDEF)
"""
abstract type AbstractNoTimeSolution{T,N} <: AbstractArray{T,N} end

"""
$(TYPEDEF)
"""
abstract type AbstractTimeseriesSolution{T,N} <: AbstractDiffEqArray{T,N} end

"""
$(TYPEDEF)
"""
abstract type AbstractEnsembleSolution{T,N} <: AbstractVectorOfArray{T,N} end

"""
$(TYPEDEF)
"""
abstract type AbstractNoiseProcess{T,N,isinplace} <: AbstractDiffEqArray{T,N} end

const DESolution = Union{AbstractTimeseriesSolution,
                         AbstractNoTimeSolution,
                         AbstractEnsembleSolution,
                         AbstractNoiseProcess}
export DESolution

"""
$(TYPEDEF)
"""
abstract type AbstractSteadyStateSolution{T,N} <: AbstractNoTimeSolution{T,N} end

"""
$(TYPEDEF)
"""
abstract type AbstractAnalyticalSolution{T,N} <: AbstractTimeseriesSolution{T,N} end

"""
$(TYPEDEF)
"""
abstract type AbstractODESolution{T,N} <: AbstractTimeseriesSolution{T,N} end

# Needed for plot recipes
"""
$(TYPEDEF)
"""
abstract type AbstractDDESolution{T,N} <: AbstractODESolution{T,N} end

"""
$(TYPEDEF)
"""
abstract type AbstractRODESolution{T,N} <: AbstractODESolution{T,N} end

"""
$(TYPEDEF)
"""
abstract type AbstractDAESolution{T,N} <: AbstractODESolution{T,N} end

"""
$(TYPEDEF)
"""
abstract type AbstractSensitivitySolution{T,N} end

# Misc
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
abstract type AbstractDiffEqOperator{T} end

"""
$(TYPEDEF)
"""
abstract type AbstractDiffEqLinearOperator{T} <: AbstractDiffEqOperator{T} end

"""
$(TYPEDEF)

Base for types defining differential equation functions.
"""
abstract type AbstractDiffEqFunction{iip} <: Function end

"""
$(TYPEDEF)

Base for types which define the history of a delay differential equation.
"""
abstract type AbstractHistoryFunction <: Function end

"""
$(TYPEDEF)
"""
abstract type AbstractReactionNetwork <: Function end

include("diffeqfastbc.jl")
include("destats.jl")
include("utils.jl")
include("calculate_residuals.jl")
include("solutions/steady_state_solutions.jl")
include("solutions/ode_solutions.jl")
include("solutions/rode_solutions.jl")
include("solutions/dae_solutions.jl")
include("solutions/solution_interface.jl")
include("tableaus.jl")
include("diffeqfunction.jl")
include("problems/problem_utils.jl")
include("problems/discrete_problems.jl")
include("problems/steady_state_problems.jl")
include("problems/analytical_problems.jl")
include("problems/ode_problems.jl")
include("problems/rode_problems.jl")
include("problems/sde_problems.jl")
include("problems/noise_problems.jl")
include("problems/bvp_problems.jl")
include("problems/dae_problems.jl")
include("problems/dde_problems.jl")
include("problems/problem_traits.jl")
include("ensemble/ensemble_solutions.jl")
include("ensemble/ensemble_problems.jl")
include("ensemble/basic_ensemble_solve.jl")
include("ensemble/ensemble_analysis.jl")
include("nlsolve/type.jl")
include("nlsolve/newton.jl")
include("nlsolve/functional.jl")
include("nlsolve/utils.jl")
include("operators/diffeq_operator.jl")
include("operators/common_defaults.jl")
include("operators/basic_operators.jl")
include("interpolation.jl")
include("callbacks.jl")
include("integrator_interface.jl")
include("linear_nonlinear.jl")
include("common_defaults.jl")
include("data_array.jl")
include("solve.jl")
include("internal_euler.jl")
include("tabletraits.jl")
include("alg_traits.jl")
include("remake.jl")
include("init.jl")

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

export isinplace

export solve, solve!, init, step!

export tuples, intervals, TimeChoiceIterator

export resize!,deleteat!,addat!,get_tmp_cache,
       full_cache,user_cache,u_cache,du_cache,
       rand_cache,ratenoise_cache,
       resize_non_user_cache!,deleteat_non_user_cache!,addat_non_user_cache!,
       terminate!,
       add_tstop!,add_saveat!,set_abstol!,
       set_reltol!,get_du, get_du!, get_dt,get_proposed_dt,set_proposed_dt!,
       u_modified!, savevalues!,reinit!, auto_dt_reset!, set_t!,
       set_u!, check_error, change_t_via_interpolation!, addsteps!,
       isdiscrete, reeval_internals_due_to_modification!

export DiscreteProblem

export SteadyStateProblem, SteadyStateSolution
export NoiseProblem
export ODEProblem, ODESolution
export DynamicalODEFunction, DynamicalODEProblem,
       SecondOrderODEProblem, SplitFunction, SplitODEProblem
export SplitSDEProblem
export RODEProblem, RODESolution, SDEProblem
export DAEProblem, DAESolution
export DDEProblem

export BVProblem, TwoPointBVProblem

export remake

export DEDataArray, DEDataVector, DEDataMatrix

export ODEFunction, DiscreteFunction, SplitFunction, DAEFunction, DDEFunction,
       SDEFunction, SplitSDEFunction, RODEFunction

export ContinuousCallback, DiscreteCallback, CallbackSet, VectorContinuousCallback

export initialize!

export LinSolveFactorize, LinSolveGPUFactorize, DefaultLinSolve, DEFAULT_LINSOLVE,
       LinSolveGMRES, LinSolveCG, LinSolveBiCGStabl, LinSolveChebyshev,
       LinSolveMINRES, LinSolveIterativeSolvers

export AffineDiffEqOperator, update_coefficients!, update_coefficients, isconstant,
       has_expmv!, has_expmv, has_exp, has_mul, has_mul!, has_ldiv, has_ldiv!

export DiffEqScalar, DiffEqArrayOperator, DiffEqIdentity

export NLNewton, NLFunctional, NLAnderson

export EnsembleThreads, EnsembleDistributed, EnsembleSplitThreads, EnsembleSerial

export EnsembleAnalysis

end # module
