__precompile__()

module DiffEqBase

  using RecipesBase, Parameters, RecursiveArrayTools
  using Ranges # For plot recipes with units
  import Base: length, size, getindex, endof, show, print,
               next, start, done, eltype, eachindex,
               Callable

  "`DEProblem`: Defines differential equation problems via its internal functions"
  abstract DEProblem
  abstract DEElement
  abstract DESensitivity
  abstract AbstractODEProblem{uType,tType,isinplace,F} <: DEProblem
  abstract AbstractODETestProblem{uType,tType,isinplace,F} <: AbstractODEProblem{uType,tType,isinplace,F}
  abstract AbstractSDEProblem <: DEProblem
  abstract AbstractDAEProblem <: DEProblem
  abstract AbstractDDEProblem <: DEProblem
  abstract AbstractPoissonProblem <: DEProblem
  abstract AbstractHeatProblem <: DEProblem
  "`PdeSolution`: Wrapper for the objects obtained from a solver"
  abstract DESolution
  abstract AbstractODESolution <: DESolution
  abstract AbstractSDESolution <: AbstractODESolution # Needed for plot recipes
  abstract AbstractDDESolution <: AbstractODESolution # Needed for plot recipes
  abstract AbstractFEMSolution <: DESolution
  abstract AbstractSensitivitySolution
  "`Mesh`: An abstract type which holds a (node,elem) pair and other information for a mesh"
  abstract Mesh
  "`Tableau`: Holds the information for a Runge-Kutta Tableau"
  abstract Tableau

  "`ODERKTableau`: A Runge-Kutta Tableau for an ODE integrator"
  abstract ODERKTableau <: Tableau

  "`DEIntegrator`: A DifferentialEquations Integrator type, used to initiate a solver."
  abstract DEIntegrator

  abstract ParameterizedFunction <: Function
  abstract SensitivityFunction <: ParameterizedFunction

  include("utils.jl")
  include("solutions/ode_solutions.jl")
  include("algorithms/ode_algorithms.jl")
  include("algorithms/ode_default_alg.jl")
  include("solutions/solution_interface.jl")
  include("tableaus.jl")
  include("problems/ode_problems.jl")
  include("problems/functions.jl")

  function solve end

  function solve(prob::DEProblem,args...;default_set=false,kwargs...)
    if default_set == true && !isempty(args)
      error("The chosen algorithm, "*string(args[1])*", does not exist.
        Please verify that the appropriate solver package has been installed.")
    elseif default_set == true
      error("No algorithm is chosen but default_set=true.")
    end
    alg,extra_kwargs = default_algorithm(prob;kwargs...)
    solve(prob,alg,args...;default_set=true,kwargs...,extra_kwargs...)
  end

  export DEProblem, DESolution, DEParameters, AbstractDAEProblem, AbstractDDEProblem,
         AbstractODEProblem, AbstractSDEProblem, DAESolution, DEIntegrator, Mesh,
         Tableau, DESensitivity, AbstractODESolution, AbstractPoissonProblem,
         AbstractHeatProblem, AbstractFEMSolution, ODERKTableau, ExplicitRKTableau,
         ImplicitRKTableau, AbstractSDESolution, ParameterizedFunction,
         AbstractSensitivitySolution, SensitivityFunction, DESensitivity, solve

  export tuples

  export @def, numparameters

  export jac_exists, invjac_exists, hes_exists, invhes_exists,
        paramjac_exists, pfunc_exists, pderiv_exists

  export ODEProblem, ODETestProblem, ODESolution, ODETestSolution,
         build_ode_solution

  # ODE Algorithms

  export OrdinaryDiffEqAlgorithm, OrdinaryDiffEqAdaptiveAlgorithm,
        Euler, Midpoint, RK4, ExplicitRK, BS3, BS5, DP5, DP5Threaded, Tsit5,
        DP8, Vern6, Vern7, Vern8, TanYam7, TsitPap8, Vern9, ImplicitEuler,
        Trapezoid, Rosenbrock23, Rosenbrock32, Feagin10, Feagin12, Feagin14

  export default_algorithm

  export ODEInterfaceAlgorithm, dopri5, dop853, odex, seulex, radau, radau5

  export ODEIterAlgorithm, feuler, rk23, feh45, feh78, ModifiedRosenbrockIntegrator,
        midpoint, heun, rk4, rk45

  export ODEJLAlgorithm, ode1, ode23, ode45, ode78, ode23s, ode2_midpoint, ode2_heun,
        ode4, ode45_fe

  export SundialsAlgorithm, CVODE_BDF, CVODE_Adams

end # module
