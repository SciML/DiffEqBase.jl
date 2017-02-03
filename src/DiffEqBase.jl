__precompile__()

module DiffEqBase

  using RecipesBase, SimpleTraits, RecursiveArrayTools
  using Ranges # For plot recipes with units
  import Base: length, size, getindex, endof, show, print,
               next, start, done, eltype, eachindex

  import Base: resize!, deleteat!

  # Problems
  "`DEProblem`: Defines differential equation problems via its internal functions"
  abstract DEProblem
  abstract DEElement
  abstract DESensitivity
  abstract AbstractODEProblem{uType,tType,isinplace,F} <: DEProblem
  abstract AbstractODETestProblem{uType,tType,isinplace,F} <: AbstractODEProblem{uType,tType,isinplace,F}
  abstract AbstractSDEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3}  <: DEProblem
  abstract AbstractSDETestProblem{uType,tType,isinplace,NoiseClass,F,F2,F3}  <: AbstractSDEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3}
  abstract AbstractDAEProblem{uType,duType,tType,isinplace,F} <: DEProblem
  abstract AbstractDAETestProblem{uType,duType,tType,isinplace,F} <: AbstractDAEProblem{uType,duType,tType,isinplace,F}
  abstract AbstractDDEProblem{uType,tType,lType,isinplace,F,H} <: DEProblem
  abstract AbstractDDETestProblem{uType,tType,lType,isinplace,F,H} <: AbstractDDEProblem{uType,tType,lType,isinplace,F,H}
  abstract AbstractConstantLagDDEProblem{uType,tType,lType,isinplace,F,H} <: AbstractDDEProblem{uType,tType,lType,isinplace,F,H}
  abstract AbstractConstantLagDDETestProblem{uType,tType,lType,isinplace,F,H} <: AbstractDDETestProblem{uType,tType,lType,isinplace,F,H}

  # Algorithms
  abstract DEAlgorithm
  abstract AbstractODEAlgorithm <: DEAlgorithm
  abstract AbstractSDEAlgorithm <: DEAlgorithm
  abstract AbstractDAEAlgorithm <: DEAlgorithm
  abstract AbstractDDEAlgorithm <: DEAlgorithm

  # Monte Carlo Simulations
  abstract AbstractMonteCarloSimulation

  # Options
  abstract DEOptions

  # Caches
  abstract DECache

  # Callbacks
  abstract DECallback

  # Integrators
  abstract DEIntegrator
  abstract AbstractODEIntegrator <: DEIntegrator
  abstract AbstractSDEIntegrator <: DEIntegrator
  abstract AbstractDDEIntegrator <: DEIntegrator
  abstract AbstractDAEIntegrator <: DEIntegrator

  # Solutions
  abstract DESolution
  abstract DETestSolution <: DESolution
  abstract AbstractODESolution <: DESolution
  abstract AbstractODETestSolution <: AbstractODESolution # Needed for plot recipes
  abstract AbstractSDESolution <: AbstractODESolution # Needed for plot recipes
  abstract AbstractSDETestSolution <: AbstractODESolution # Needed for plot recipes
  abstract AbstractDAESolution <: AbstractODESolution # Needed for plot recipes
  abstract AbstractDAETestSolution <: AbstractODESolution # Needed for plot recipes
  abstract AbstractDDESolution <: DESolution
  abstract AbstractDDETestSolution <: AbstractODESolution # Needed for plot recipes
  abstract AbstractSensitivitySolution

  # Misc
  "`Tableau`: Holds the information for a Runge-Kutta Tableau"
  abstract Tableau
  "`ODERKTableau`: A Runge-Kutta Tableau for an ODE integrator"
  abstract ODERKTableau <: Tableau

  abstract AbstractParameterizedFunction{isinplace} <: Function

  include("utils.jl")
  include("extended_functions.jl")
  include("noise_process.jl")
  include("solutions/ode_solutions.jl")
  include("solutions/sde_solutions.jl")
  include("solutions/dae_solutions.jl")
  include("solutions/dde_solutions.jl")
  include("solutions/solution_interface.jl")
  include("tableaus.jl")
  include("problems/ode_problems.jl")
  include("problems/sde_problems.jl")
  include("problems/dae_problems.jl")
  include("problems/dde_problems.jl")
  include("callbacks.jl")
  include("integrator_interface.jl")

  type ConvergenceSetup{P,C}
    probs::P
    convergence_axis::C
  end

  function solve end
  function solve! end
  function init end
  function step!(d::DEIntegrator) error("Integrator stepping is not implemented") end

  export DEProblem, DESolution, DETestSolution, DEParameters, AbstractDAEProblem,
         AbstractDDEProblem, AbstractODEProblem, AbstractSDEProblem, DAESolution,
         DEIntegrator, Mesh, Tableau, DESensitivity, AbstractODESolution, ODERKTableau,
         ExplicitRKTableau,
         ImplicitRKTableau, AbstractSDESolution, DESensitivity, DEAlgorithm,
         AbstractODETestProblem, DECallback, DECache, DEIntegrator,
         DEOptions, AbstractODETestSolution, AbstractDDEProblem, AbstractDDETestProblem,
         AbstractSDETestProblem, AbstractDAETestProblem,   AbstractSDETestSolution,
         AbstractDAETestSolution, AbstractDDESolution, AbstractDDETestSolution,
         AbstractMonteCarloSimulation, AbstractConstantLagDDETestProblem,
         AbstractConstantLagDDEProblem

  export isinplace

  export solve, solve!, init, step!

  export tuples, intervals

  export resize!,deleteat!,full_cache,u_cache,du_cache,terminate!,add_tstop!,add_saveat!,set_abstol!,
         set_reltol!,get_du,get_dt,get_proposed_dt,modify_proposed_dt!,u_modified!,
         savevalues!

  export numparameters, @def, @muladd

  export NoiseProcess, construct_correlated_noisefunc

  export HasJac, HastGrad, HasParamFuncs, HasParamDeriv, HasParamJac,
         HasInvJac,HasInvW, HasInvW_t, HasHes, HasInvHes, HasSyms

  export has_jac, has_invjac, has_invW, has_invW_t, has_hes, has_invhes,
         has_tgrad, has_paramfuncs, has_paramderiv, has_paramjac,
         has_syms

  export ODEProblem, ODETestProblem, ODESolution, ODETestSolution

  export SDEProblem, SDETestProblem, SDESolution, SDETestSolution

  export DAEProblem, DAETestProblem, DAESolution, DAETestSolution

  export ConstantLagDDEProblem, ConstantLagDDETestProblem, DDEProblem, DDETestProblem

  export build_solution, calculate_solution_errors!

  export AbstractParameterizedFunction

  export ConvergenceSetup

  export ContinuousCallback, DiscreteCallback, CallbackSet

  # Algorithms

  export AbstractODEAlgorithm, AbstractSDEAlgorithm, AbstractDAEAlgorithm, AbstractDDEAlgorithm

  export AbstractODEIntegrator, AbstractSDEIntegrator, AbstractDDEIntegrator, AbstractDAEIntegrator


end # module
