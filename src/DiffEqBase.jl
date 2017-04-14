__precompile__()

module DiffEqBase

using RecipesBase, SimpleTraits, RecursiveArrayTools, Compat
using Ranges # For plot recipes with units
import Base: length, ndims, size, getindex, setindex!, endof, show, print,
             next, start, done, eltype, eachindex, similar

import Base: resize!, deleteat!

import RecursiveArrayTools: recursivecopy!

# Problems
@compat abstract type DEProblem end
@compat abstract type DEElement end
@compat abstract type DESensitivity end

export DEProblem, DEElement, DESensitivity


@compat abstract type AbstractODEProblem{uType,tType,isinplace} <: DEProblem end
@compat abstract type AbstractDiscreteProblem{uType,tType,isinplace} <:
                      AbstractODEProblem{uType,tType,isinplace} end

export AbstractODEProblem, AbstractDiscreteProblem

@compat abstract type AbstractRODEProblem{uType,tType,isinplace,NoiseClass,ND} <: DEProblem end

export AbstractRODEProblem

@compat abstract type AbstractSDEProblem{uType,tType,isinplace,NoiseClass,ND} <:
                      AbstractRODEProblem{uType,tType,isinplace,NoiseClass,ND} end

export AbstractSDEProblem

@compat abstract type AbstractDAEProblem{uType,duType,tType,isinplace} <: DEProblem end

export AbstractDAEProblem

@compat abstract type AbstractDDEProblem{uType,tType,lType,isinplace} <: DEProblem end
@compat abstract type AbstractConstantLagDDEProblem{uType,tType,lType,isinplace} <:
                      AbstractDDEProblem{uType,tType,lType,isinplace} end

export AbstractDDEProblem,
       AbstractConstantLagDDEProblem

@compat abstract type AbstractSecondOrderODEProblem{uType,tType,isinplace} <: AbstractODEProblem{uType,tType,isinplace} end

export AbstractSecondOrderODEProblem

# Algorithms
@compat abstract type DEAlgorithm end
@compat abstract type AbstractODEAlgorithm <: DEAlgorithm end
@compat abstract type AbstractSplitODEAlgorithm <: DEAlgorithm end
@compat abstract type AbstractPartitionedODEAlgorithm <: DEAlgorithm end
@compat abstract type AbstractSecondOrderODEAlgorithm <: DEAlgorithm end
@compat abstract type AbstractSplitDAEAlgorithm <: DEAlgorithm end
@compat abstract type AbstractPartitionedDAEAlgorithm <: DEAlgorithm end
@compat abstract type AbstractSplitConstrainedODEAlgorithm <: DEAlgorithm end
@compat abstract type AbstractPartitionedConstrainedODEAlgorithm <: DEAlgorithm end
@compat abstract type AbstractRODEAlgorithm <: DEAlgorithm end
@compat abstract type AbstractSDEAlgorithm <: DEAlgorithm end
@compat abstract type AbstractDAEAlgorithm <: DEAlgorithm end
@compat abstract type AbstractDDEAlgorithm <: DEAlgorithm end

export DEAlgorithm, AbstractODEAlgorithm, AbstractSplitODEAlgorithm,
       AbstractPartitionedODEAlgorithm, AbstractSecondOrderODEAlgorithm,
       AbstractRODEAlgorithm, AbstractSDEAlgorithm,
       AbstractDAEAlgorithm, AbstractDDEAlgorithm

# Monte Carlo Simulations
@compat abstract type AbstractMonteCarloProblem <: DEProblem end
@compat abstract type AbstractMonteCarloEstimator <: DEProblem end

export AbstractMonteCarloProblem, MonteCarloProblem
export MonteCarloSolution, MonteCarloTestSolution
export AbstractMonteCarloEstimator

@compat abstract type DEOptions end
@compat abstract type DECache end
@compat abstract type DECallback end
@compat abstract type DEDataArray{T} <: AbstractVector{T} end

export DEOptions, DECache, DECallback, DEDataArray

# Integrators
@compat abstract type DEIntegrator end
@compat abstract type AbstractODEIntegrator <: DEIntegrator end
@compat abstract type AbstractSplitODEIntegrator <: DEIntegrator end
@compat abstract type AbstractPartitionedODEIntegrator <: DEIntegrator end
@compat abstract type AbstractSecondOrderODEIntegrator <: DEIntegrator end
@compat abstract type AbstractRODEIntegrator <: DEIntegrator end
@compat abstract type AbstractSDEIntegrator <: DEIntegrator end
@compat abstract type AbstractDDEIntegrator <: DEIntegrator end
@compat abstract type AbstractDAEIntegrator <: DEIntegrator end

export DEIntegrator, AbstractODEIntegrator, AbstractSplitODEIntegrator,
       AbstractPartitionedODEIntegrator, AbstractSecondOrderODEIntegrator,
       AbstractRODEIntegrator, AbstractSDEIntegrator,
       AbstractDDEIntegrator, AbstractDAEIntegrator

# Solutions
@compat abstract type DESolution end
@compat abstract type AbstractODESolution <: DESolution end

export DESolution, AbstractODESolution, AbstractDDESolution

# Needed for plot recipes
@compat abstract type AbstractDDESolution <: AbstractODESolution end
@compat abstract type  AbstractRODESolution <: AbstractODESolution end
@compat abstract type  AbstractDAESolution <: AbstractODESolution end
@compat abstract type AbstractSensitivitySolution end
@compat abstract type AbstractMonteCarloSolution <: DESolution end

export AbstractRODESolution,
       AbstractDAESolution,
       AbstractDDESolution,
       AbstractSensitivitySolution, AbstractMonteCarloSolution

# Misc
@compat abstract type Tableau end
@compat abstract type ODERKTableau <: Tableau end
@compat abstract type DECostFunction end

export Tableau, ODERKTableau, DECostFunction

@compat abstract type AbstractParameterizedFunction{isinplace} <: Function end
@compat abstract type ConstructedParameterizedFunction{isinplace} <: AbstractParameterizedFunction{isinplace} end
export AbstractParameterizedFunction, ConstructedParameterizedFunction

include("utils.jl")
include("extended_functions.jl")
include("noise_process.jl")
include("solutions/ode_solutions.jl")
include("solutions/rode_solutions.jl")
include("solutions/dae_solutions.jl")
include("solutions/monte_solutions.jl")
include("solutions/solution_interface.jl")
include("tableaus.jl")
include("problems/discrete_problems.jl")
include("problems/ode_problems.jl")
include("problems/constrained_ode_problems.jl")
include("problems/rode_problems.jl")
include("problems/sde_problems.jl")
include("problems/dae_problems.jl")
include("problems/dde_problems.jl")
include("problems/monte_problems.jl")
include("problems/problem_traits.jl")
include("callbacks.jl")
include("integrator_interface.jl")
include("parameters_interface.jl")
include("constructed_parameterized_functions.jl")
include("data_array.jl")

type ConvergenceSetup{P,C}
    probs::P
    convergence_axis::C
end

function solve end
function solve! end
function init end
function step!(d::DEIntegrator) error("Integrator stepping is not implemented") end

export DEParameters, Mesh, ExplicitRKTableau, ImplicitRKTableau

export isinplace, noise_class

export solve, solve!, init, step!

export tuples, intervals

export resize!,deleteat!,addat!,full_cache,user_cache,u_cache,du_cache,
       resize_non_user_cache!,deleteat_non_user_cache!,addat_non_user_cache!,
       terminate!,
       add_tstop!,add_saveat!,set_abstol!,
       set_reltol!,get_du,get_dt,get_proposed_dt,set_proposed_dt!,u_modified!,
       savevalues!

export numargs, @def, @muladd

export NoiseProcess, construct_correlated_noisefunc

export HasJac, HastGrad, HasParamFuncs, HasParamDeriv, HasParamJac,
       HasInvJac,HasInvW, HasInvW_t, HasHes, HasInvHes, HasSyms

export has_jac, has_invjac, has_invW, has_invW_t, has_hes, has_invhes,
       has_tgrad, has_paramfuncs, has_paramderiv, has_paramjac,
       has_syms, has_analytic

export DiscreteProblem

export ODEProblem, ODESolution
export SplitODEProblem, PartitionedODEProblem, SecondOrderODEProblem
export RODEProblem, RODESolution, SDEProblem
export DAEProblem, DAESolution
export SplitDAEProblem, PartitionedDAEProblem, ConstrainedODEProblem,
       SplitConstrainedODEProblem, PartitionedConstrainedODEProblem
export ConstantLagDDEProblem, DDEProblem

export ParameterizedFunction, DAEParameterizedFunction, DDEParameterizedFunction,
       ProbParameterizedFunction, OutputParameterizedFunction

export calculate_monte_errors

export build_solution, calculate_solution_errors!

export ConvergenceSetup

export ContinuousCallback, DiscreteCallback, CallbackSet

export initialize!

export copy_non_array_fields

export param_values, num_params, problem_new_parameters

export white_noise_func_wrapper, white_noise_func_wrapper!

export is_diagonal_noise, is_sparse_noise

end # module
