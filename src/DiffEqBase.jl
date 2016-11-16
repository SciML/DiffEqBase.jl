__precompile__()

module DiffEqBase

  using RecipesBase, RecursiveArrayTools
  using Ranges # For plot recipes with units
  import Base: length, size, getindex, endof, show, print,
               next, start, done, eltype, eachindex

  "`DEProblem`: Defines differential equation problems via its internal functions"
  abstract DEProblem
  abstract DEElement
  abstract DESensitivity
  abstract AbstractODEProblem{uType,tType,isinplace,F} <: DEProblem
  abstract AbstractODETestProblem{uType,tType,isinplace,F} <: AbstractODEProblem{uType,tType,isinplace,F}
  abstract AbstractSDEProblem <: DEProblem
  abstract AbstractDAEProblem <: DEProblem
  abstract AbstractDDEProblem <: DEProblem
  abstract AbstractPDEProblem <: DEProblem
  "`DSolution`: type to hold the objects obtained from a solver"
  abstract DESolution
  abstract AbstractODESolution <: DESolution
  abstract AbstractSDESolution <: AbstractODESolution # Needed for plot recipes
  abstract AbstractDDESolution <: AbstractODESolution # Needed for plot recipes
  abstract AbstractFEMSolution <: DESolution
  abstract AbstractSensitivitySolution

"`Tableau`: Holds the information for a Runge-Kutta Tableau"
  abstract Tableau

  "`ODERKTableau`: A Runge-Kutta Tableau for an ODE integrator"
  abstract ODERKTableau <: Tableau

  include("utils.jl")
  include("solutions/ode_solutions.jl")
  include("solutions/solution_interface.jl")
  include("tableaus.jl")
  include("problems/ode_problems.jl")

  function solve end

  # function solve(prob::DEProblem,args...;default_set=false,kwargs...)
  #   if default_set == true && !isempty(args)
  #     error("The chosen algorithm, "*string(args[1])*", does not exist.
  #       Please verify that the appropriate solver package has been installed.")
  #   elseif default_set == true
  #     error("No algorithm is chosen but default_set=true.")
  #   end
  #   alg,extra_kwargs = default_algorithm(prob;kwargs...)
  #   solve(prob,alg,args...;default_set=true,kwargs...,extra_kwargs...)
  # end

  export DEProblem, DESolution, DEParameters, AbstractDAEProblem, AbstractDDEProblem,
         AbstractODEProblem, AbstractSDEProblem, DAESolution,
         Tableau, DESensitivity, AbstractODESolution, AbstractPoissonProblem,
         AbstractHeatProblem, AbstractFEMSolution, ODERKTableau, ExplicitRKTableau,
         ImplicitRKTableau, AbstractSDESolution, ParameterizedFunction,
         AbstractSensitivitySolution, SensitivityFunction, DESensitivity, solve

  export tuples

  export numparameters

  export jac_exists, invjac_exists, hes_exists, invhes_exists,
        paramjac_exists, pfunc_exists, pderiv_exists

  export ODEProblem, ODETestProblem, ODESolution, ODETestSolution,
         build_ode_solution

end # module
