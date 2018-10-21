module DiffEqBase

using RecipesBase, RecursiveArrayTools, Compat,
      Requires, TableTraits, IteratorInterfaceExtensions, TreeViews

import Base: length, size, getindex, setindex!, show, print,
             iterate, eltype, similar, axes

import Base: resize!, deleteat!

import RecursiveArrayTools: recursivecopy!, tuples

import RecursiveArrayTools: chain

import Roots # callbacks

import LinearAlgebra: exp, ldiv!

import StaticArrays: StaticArray, SArray, Size

import Statistics: mean, median

# Problems
abstract type DEProblem end
abstract type DEElement end
abstract type DESensitivity end

abstract type AbstractSteadyStateProblem{uType,isinplace} <: DEProblem end

abstract type AbstractNoiseProblem <: DEProblem end

abstract type AbstractODEProblem{uType,tType,isinplace} <: DEProblem end
abstract type AbstractDiscreteProblem{uType,tType,isinplace} <:
                      AbstractODEProblem{uType,tType,isinplace} end
abstract type AbstractAnalyticalProblem{uType,tType,isinplace} <:
                      AbstractODEProblem{uType,tType,isinplace} end

abstract type AbstractRODEProblem{uType,tType,isinplace,ND} <: DEProblem end

abstract type AbstractSDEProblem{uType,tType,isinplace,ND} <:
                      AbstractRODEProblem{uType,tType,isinplace,ND} end

abstract type AbstractDAEProblem{uType,duType,tType,isinplace} <: DEProblem end

abstract type AbstractDDEProblem{uType,tType,lType,isinplace} <: DEProblem end
abstract type AbstractConstantLagDDEProblem{uType,tType,lType,isinplace} <:
                      AbstractDDEProblem{uType,tType,lType,isinplace} end

abstract type AbstractSecondOrderODEProblem{uType,tType,isinplace} <: AbstractODEProblem{uType,tType,isinplace} end

abstract type AbstractBVProblem{uType,tType,isinplace} <: AbstractODEProblem{uType,tType,isinplace} end

abstract type AbstractJumpProblem{P,J} <: DiffEqBase.DEProblem end

# Algorithms
abstract type DEAlgorithm end
abstract type AbstractSteadyStateAlgorithm <: DEAlgorithm end
abstract type AbstractODEAlgorithm <: DEAlgorithm end
abstract type AbstractSecondOrderODEAlgorithm <: DEAlgorithm end
abstract type AbstractRODEAlgorithm <: DEAlgorithm end
abstract type AbstractSDEAlgorithm <: DEAlgorithm end
abstract type AbstractDAEAlgorithm <: DEAlgorithm end
abstract type AbstractDDEAlgorithm <: DEAlgorithm end

# Monte Carlo Simulations
abstract type AbstractMonteCarloProblem <: DEProblem end
abstract type AbstractMonteCarloEstimator <: DEProblem end

export MonteCarloProblem
export MonteCarloSolution, MonteCarloTestSolution, MonteCarloSummary

abstract type AbstractDiffEqInterpolation <: Function end
abstract type AbstractDEOptions end
abstract type DECache end
abstract type DECallback end
abstract type AbstractContinuousCallback <: DECallback end
abstract type AbstractDiscreteCallback <: DECallback end

abstract type DEDataArray{T,N} <: AbstractArray{T,N} end
const DEDataVector{T} = DEDataArray{T,1}
const DEDataMatrix{T} = DEDataArray{T,2}

# Integrators
abstract type DEIntegrator end
abstract type AbstractSteadyStateIntegrator <: DEIntegrator end
abstract type AbstractODEIntegrator <: DEIntegrator end
abstract type AbstractSecondOrderODEIntegrator <: DEIntegrator end
abstract type AbstractRODEIntegrator <: DEIntegrator end
abstract type AbstractSDEIntegrator <: DEIntegrator end
abstract type AbstractDDEIntegrator <: DEIntegrator end
abstract type AbstractDAEIntegrator <: DEIntegrator end

# Solutions
abstract type AbstractNoTimeSolution{T,N} <: AbstractArray{T,N} end
abstract type AbstractTimeseriesSolution{T,N} <: AbstractDiffEqArray{T,N} end
abstract type AbstractMonteCarloSolution{T,N} <: AbstractVectorOfArray{T,N} end
abstract type AbstractNoiseProcess{T,N,isinplace} <: AbstractDiffEqArray{T,N} end

const DESolution = Union{AbstractTimeseriesSolution,
                         AbstractNoTimeSolution,
                         AbstractMonteCarloSolution,
                         AbstractNoiseProcess}
export DESolution
abstract type AbstractSteadyStateSolution{T,N} <: AbstractNoTimeSolution{T,N} end
abstract type AbstractAnalyticalSolution{T,N} <: AbstractTimeseriesSolution{T,N} end
abstract type AbstractODESolution{T,N} <: AbstractTimeseriesSolution{T,N} end

# Needed for plot recipes
abstract type AbstractDDESolution{T,N} <: AbstractODESolution{T,N} end
abstract type AbstractRODESolution{T,N} <: AbstractODESolution{T,N} end
abstract type AbstractDAESolution{T,N} <: AbstractODESolution{T,N} end
abstract type AbstractSensitivitySolution{T,N} end

# Misc
abstract type Tableau end
abstract type ODERKTableau <: Tableau end
abstract type DECostFunction end
abstract type AbstractDiffEqOperator{T} end
abstract type AbstractDiffEqLinearOperator{T} <: AbstractDiffEqOperator{T} end
abstract type AbstractDiffEqFunction{iip} <: Function end
abstract type AbstractReactionNetwork <: Function end

include("utils.jl")
include("extended_functions.jl")
include("solutions/steady_state_solutions.jl")
include("solutions/ode_solutions.jl")
include("solutions/rode_solutions.jl")
include("solutions/dae_solutions.jl")
include("solutions/monte_solutions.jl")
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
include("problems/monte_problems.jl")
include("problems/problem_traits.jl")
include("interpolation.jl")
include("callbacks.jl")
include("integrator_interface.jl")
include("parameters_interface.jl")
include("diffeq_operator.jl")
include("linear_nonlinear.jl")
include("common_defaults.jl")
include("data_array.jl")
include("solve.jl")
include("internal_euler.jl")
include("tabletraits.jl")
include("alg_traits.jl")
include("remake.jl")
include("init.jl")

abstract type AbstractParameterizedFunction{iip} <: AbstractODEFunction{iip} end

struct ConvergenceSetup{P,C}
    probs::P
    convergence_axis::C
end

export isinplace

export solve, solve!, init, step!

export tuples, intervals, TimeChoiceIterator

export resize!,deleteat!,addat!,get_tmp_cache,full_cache,user_cache,u_cache,du_cache,
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

export ContinuousCallback, DiscreteCallback, CallbackSet

export initialize!

export problem_new_parameters

export LinSolveFactorize, DEFAULT_LINSOLVE

export AffineDiffEqOperator, update_coefficients!, update_coefficients, is_constant,
       has_expmv!, has_expmv, has_exp, has_mul, has_mul!, has_ldiv, has_ldiv!

end # module
