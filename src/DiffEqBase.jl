__precompile__()

module DiffEqBase

  using RecipesBase, Parameters, RecursiveArrayTools
  using Ranges # For plot recipes with units
  import Base: length, size, getindex, endof, show, print,
               next, start, done, eltype, eachindex,
               Callable

  # Problems
  "`DEProblem`: Defines differential equation problems via its internal functions"
  abstract DEProblem
  abstract DEElement
  abstract DESensitivity
  abstract AbstractODEProblem{uType,tType,isinplace,F} <: DEProblem
  abstract AbstractODETestProblem{uType,tType,isinplace,F} <: AbstractODEProblem{uType,tType,isinplace,F}
  abstract AbstractSDEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3}  <: DEProblem
  abstract AbstractSDETestProblem{uType,tType,isinplace,NoiseClass,F,F2,F3}  <: AbstractSDEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3}
  abstract AbstractDAEProblem <: DEProblem
  abstract AbstractDDEProblem <: DEProblem
  abstract AbstractPoissonProblem <: DEProblem
  abstract AbstractHeatProblem <: DEProblem

  # Algorithms
  abstract DEAlgorithm
  abstract AbstractODEAlgorithm <: DEAlgorithm
  abstract AbstractSDEAlgorithm <: DEAlgorithm

  # Solutions
  abstract DESolution
  abstract AbstractODESolution <: DESolution
  abstract AbstractSDESolution <: AbstractODESolution # Needed for plot recipes
  abstract AbstractDDESolution <: AbstractODESolution # Needed for plot recipes
  abstract AbstractFEMSolution <: DESolution
  abstract AbstractSensitivitySolution

  # Misc
  "`Mesh`: An abstract type which holds a (node,elem) pair and other information for a mesh"
  abstract Mesh
  "`Tableau`: Holds the information for a Runge-Kutta Tableau"
  abstract Tableau
  "`ODERKTableau`: A Runge-Kutta Tableau for an ODE integrator"
  abstract ODERKTableau <: Tableau

  include("utils.jl")
  include("noise_process.jl")
  include("solutions/ode_solutions.jl")
  include("solutions/sde_solutions.jl")
  include("solutions/solution_interface.jl")
  include("tableaus.jl")
  include("problems/ode_problems.jl")
  include("problems/sde_problems.jl")
  include("problems/functions.jl")

  function solve end

  export DEProblem, DESolution, DEParameters, AbstractDAEProblem, AbstractDDEProblem,
         AbstractODEProblem, AbstractSDEProblem, DAESolution, DEIntegrator, Mesh,
         Tableau, DESensitivity, AbstractODESolution, AbstractPoissonProblem,
         AbstractHeatProblem, AbstractFEMSolution, ODERKTableau, ExplicitRKTableau,
         ImplicitRKTableau, AbstractSDESolution, DESensitivity, solve

  export tuples

  export numparameters

  export jac_exists, invjac_exists, hes_exists, invhes_exists,
        paramjac_exists, pfunc_exists, pderiv_exists

  export ODEProblem, ODETestProblem, ODESolution, ODETestSolution,
         build_ode_solution

  export SDEProblem, SDETestProblem, SDESolution, SDETestSolution,
         build_sde_solution

  # Algorithms

  export AbstractODEAlgorithm, AbstractSDEAlgorithm

end # module
