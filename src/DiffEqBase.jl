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


@compat abstract type AbstractODEProblem{uType,tType,isinplace,F} <: DEProblem end
@compat abstract type AbstractODETestProblem{uType,tType,isinplace,F} <:
                      AbstractODEProblem{uType,tType,isinplace,F} end
@compat abstract type AbstractDiscreteProblem{uType,tType,isinplace,F} <:
                      AbstractODEProblem{uType,tType,isinplace,F} end
@compat abstract type AbstractDiscreteTestProblem{uType,tType,isinplace,F} <:
                      AbstractODETestProblem{uType,tType,isinplace,F} end

export AbstractODEProblem, AbstractODETestProblem, AbstractDiscreteProblem, AbstractDiscreteTestProblem

@compat abstract type AbstractSDEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3} <: DEProblem end
@compat abstract type AbstractSDETestProblem{uType,tType,isinplace,NoiseClass,F,F2,F3} <:
                      AbstractSDEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3} end

export AbstractSDEProblem, AbstractSDETestProblem

@compat abstract type AbstractDAEProblem{uType,duType,tType,isinplace,F} <: DEProblem end
@compat abstract type AbstractDAETestProblem{uType,duType,tType,isinplace,F} <:
                      AbstractDAEProblem{uType,duType,tType,isinplace,F} end

export AbstractDAEProblem, AbstractDAETestProblem

@compat abstract type AbstractDDEProblem{uType,tType,lType,isinplace,F,H} <: DEProblem end
@compat abstract type AbstractDDETestProblem{uType,tType,lType,isinplace,F,H} <:
                      AbstractDDEProblem{uType,tType,lType,isinplace,F,H} end
@compat abstract type AbstractConstantLagDDEProblem{uType,tType,lType,isinplace,F,H} <:
                      AbstractDDEProblem{uType,tType,lType,isinplace,F,H} end
@compat abstract type AbstractConstantLagDDETestProblem{uType,tType,lType,isinplace,F,H} <:
                      AbstractDDETestProblem{uType,tType,lType,isinplace,F,H} end

export AbstractDDEProblem, AbstractDDETestProblem, AbstractConstantLagDDEProblem, AbstractConstantLagDDETestProblem

@compat abstract type AbstractSplitODEProblem{uType,tType,isinplace,F} <: AbstractODEProblem{uType,tType,isinplace,F} end
@compat abstract type AbstractSplitODETestProblem{uType,tType,isinplace,F} <: AbstractODETestProblem{uType,tType,isinplace,F} end

export AbstractSplitODEProblem, SplitODEProblem, SplitODETestProblem

# Algorithms
@compat abstract type DEAlgorithm end
@compat abstract type AbstractODEAlgorithm <: DEAlgorithm end
@compat abstract type AbstractSDEAlgorithm <: DEAlgorithm end
@compat abstract type AbstractDAEAlgorithm <: DEAlgorithm end
@compat abstract type AbstractDDEAlgorithm <: DEAlgorithm end

export DEAlgorithm, AbstractODEAlgorithm, AbstractSDEAlgorithm,
       AbstractDAEAlgorithm, AbstractDDEAlgorithm

# Monte Carlo Simulations
@compat abstract type AbstractMonteCarloProblem <: DEProblem end

export AbstractMonteCarloProblem, MonteCarloProblem
export MonteCarloSolution, MonteCarloTestSolution

@compat abstract type DEOptions end
@compat abstract type DECache end
@compat abstract type DECallback end
@compat abstract type DEDataArray{T} <: AbstractVector{T} end

export DEOptions, DECache, DECallback, DEDataArray

# Integrators
@compat abstract type DEIntegrator end
@compat abstract type AbstractODEIntegrator <: DEIntegrator end
@compat abstract type AbstractSDEIntegrator <: DEIntegrator end
@compat abstract type AbstractDDEIntegrator <: DEIntegrator end
@compat abstract type AbstractDAEIntegrator <: DEIntegrator end

export DEIntegrator, AbstractODEIntegrator, AbstractSDEIntegrator, AbstractDDEIntegrator,
       AbstractDAEIntegrator

# Solutions
@compat abstract type DESolution end
@compat abstract type DETestSolution <: DESolution end
@compat abstract type AbstractODESolution <: DESolution end

export DESolution, DETestSolution, AbstractODESolution, AbstractDDESolution

# Needed for plot recipes
@compat abstract type AbstractDDESolution <: AbstractODESolution end
@compat abstract type  AbstractODETestSolution <: AbstractODESolution end
@compat abstract type  AbstractSDESolution <: AbstractODESolution end
@compat abstract type  AbstractSDETestSolution <: AbstractODESolution end
@compat abstract type  AbstractDAESolution <: AbstractODESolution end
@compat abstract type  AbstractDAETestSolution <: AbstractODESolution end
@compat abstract type  AbstractDDETestSolution <: AbstractODESolution end
@compat abstract type AbstractSensitivitySolution end
@compat abstract type AbstractMonteCarloSolution <: DESolution end

export AbstractODETestSolution, AbstractSDESolution, AbstractSDETestSolution,
       AbstractDAESolution, AbstractDAETestSolution, AbstractDDETestSolution,
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
include("solutions/sde_solutions.jl")
include("solutions/dae_solutions.jl")
include("solutions/dde_solutions.jl")
include("solutions/monte_solutions.jl")
include("solutions/solution_interface.jl")
include("tableaus.jl")
include("problems/discrete_problems.jl")
include("problems/ode_problems.jl")
include("problems/splitode_problem.jl")
include("problems/sde_problems.jl")
include("problems/dae_problems.jl")
include("problems/dde_problems.jl")
include("problems/monte_problems.jl")
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

export isinplace

export solve, solve!, init, step!

export tuples, intervals

export resize!,deleteat!,addat!,full_cache,user_cache,u_cache,du_cache,
       resize_non_user_cache!,deleteat_non_user_cache!,addat_non_user_cache!,
       terminate!,
       add_tstop!,add_saveat!,set_abstol!,
       set_reltol!,get_du,get_dt,get_proposed_dt,set_proposed_dt!,u_modified!,
       savevalues!

export numparameters, @def, @muladd

export NoiseProcess, construct_correlated_noisefunc

export HasJac, HastGrad, HasParamFuncs, HasParamDeriv, HasParamJac,
       HasInvJac,HasInvW, HasInvW_t, HasHes, HasInvHes, HasSyms

export has_jac, has_invjac, has_invW, has_invW_t, has_hes, has_invhes,
       has_tgrad, has_paramfuncs, has_paramderiv, has_paramjac,
       has_syms

export DiscreteProblem, DiscreteTestProblem

export ODEProblem, ODETestProblem, ODESolution, ODETestSolution
export SDEProblem, SDETestProblem, SDESolution, SDETestSolution
export DAEProblem, DAETestProblem, DAESolution, DAETestSolution
export ConstantLagDDEProblem, ConstantLagDDETestProblem, DDEProblem, DDETestProblem

export ParameterizedFunction, DAEParameterizedFunction, DDEParameterizedFunction,
       ProbParameterizedFunction, OutputParameterizedFunction

export calculate_monte_errors

export build_solution, calculate_solution_errors!

export ConvergenceSetup

export ContinuousCallback, DiscreteCallback, CallbackSet

export initialize!

export copy_non_array_fields

export param_values, num_params, problem_new_parameters

end # module
