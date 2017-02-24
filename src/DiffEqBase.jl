__precompile__()

module DiffEqBase

using RecipesBase, SimpleTraits, RecursiveArrayTools
using Ranges # For plot recipes with units
import Base: length, ndims, size, getindex, setindex!, endof, show, print,
             next, start, done, eltype, eachindex, similar

import Base: resize!, deleteat!

import RecursiveArrayTools: recursivecopy!

include("typemacros.jl")
using .TypeMacros

# Problems
@public_abstract_type DEProblem
@public_abstract_type DEElement
@public_abstract_type DESensitivity


@public_abstract_type AbstractODEProblem{uType,tType,isinplace,F} <: DEProblem
@public_abstract_type AbstractODETestProblem{uType,tType,isinplace,F} <:
                      AbstractODEProblem{uType,tType,isinplace,F}
@public_abstract_type AbstractDiscreteProblem{uType,tType,isinplace,F} <:
                      AbstractODEProblem{uType,tType,isinplace,F}
@public_abstract_type AbstractDiscreteTestProblem{uType,tType,isinplace,F} <:
                      AbstractODETestProblem{uType,tType,isinplace,F}

@public_abstract_type AbstractSDEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3} <: DEProblem
@public_abstract_type AbstractSDETestProblem{uType,tType,isinplace,NoiseClass,F,F2,F3} <:
                      AbstractSDEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3}

@public_abstract_type AbstractDAEProblem{uType,duType,tType,isinplace,F} <: DEProblem
@public_abstract_type AbstractDAETestProblem{uType,duType,tType,isinplace,F} <:
                      AbstractDAEProblem{uType,duType,tType,isinplace,F}

@public_abstract_type AbstractDDEProblem{uType,tType,lType,isinplace,F,H} <: DEProblem
@public_abstract_type AbstractDDETestProblem{uType,tType,lType,isinplace,F,H} <:
                      AbstractDDEProblem{uType,tType,lType,isinplace,F,H}
@public_abstract_type AbstractConstantLagDDEProblem{uType,tType,lType,isinplace,F,H} <:
                      AbstractDDEProblem{uType,tType,lType,isinplace,F,H}
@public_abstract_type AbstractConstantLagDDETestProblem{uType,tType,lType,isinplace,F,H} <:
                      AbstractDDETestProblem{uType,tType,lType,isinplace,F,H}

# Algorithms
@public_abstract_type DEAlgorithm
@public_abstract_subs (AbstractODEAlgorithm,
                       AbstractSDEAlgorithm,
                       AbstractDAEAlgorithm,
                       AbstractDDEAlgorithm) <: DEAlgorithm

# Monte Carlo Simulations
@public_abstract_type AbstractMonteCarloSimulation

# Options
@public_abstract_type DEOptions

# Caches
@public_abstract_type DECache

# Callbacks
@public_abstract_type DECallback

# Array
@public_abstract_type DEDataArray{T} <: AbstractVector{T}

# Integrators
@public_abstract_type DEIntegrator
@public_abstract_subs (AbstractODEIntegrator,
                       AbstractSDEIntegrator,
                       AbstractDDEIntegrator,
                       AbstractDAEIntegrator) <: DEIntegrator

# Solutions
@public_abstract_type DESolution
@public_abstract_subs (DETestSolution, AbstractODESolution, AbstractDDESolution) <: DESolution

# Needed for plot recipes
@public_abstract_subs (AbstractODETestSolution,
                       AbstractSDESolution,
                       AbstractSDETestSolution,
                       AbstractDAESolution,
                       AbstractDAETestSolution,
                       AbstractDDETestSolution) <: AbstractODESolution
@abstract_type AbstractSensitivitySolution

# Misc
@public_abstract_type Tableau

@public_abstract_type ODERKTableau <: Tableau

@public_abstract_type AbstractParameterizedFunction{isinplace} <: Function

include("utils.jl")
include("extended_functions.jl")
include("noise_process.jl")
include("solutions/ode_solutions.jl")
include("solutions/sde_solutions.jl")
include("solutions/dae_solutions.jl")
include("solutions/dde_solutions.jl")
include("solutions/solution_interface.jl")
include("tableaus.jl")
include("problems/discrete_problems.jl")
include("problems/ode_problems.jl")
include("problems/sde_problems.jl")
include("problems/dae_problems.jl")
include("problems/dde_problems.jl")
include("callbacks.jl")
include("integrator_interface.jl")
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
       set_reltol!,get_du,get_dt,get_proposed_dt,modify_proposed_dt!,u_modified!,
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

export build_solution, calculate_solution_errors!

export ConvergenceSetup

export ContinuousCallback, DiscreteCallback, CallbackSet

end # module
