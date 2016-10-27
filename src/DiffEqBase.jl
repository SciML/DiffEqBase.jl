module DiffEqBase

  using RecipesBase
  using Ranges # For plot recipes with units
  import Base: length, size, getindex, endof, show, print

  "`DEProblem`: Defines differential equation problems via its internal functions"
  abstract DEProblem
  abstract DEElement
  abstract AbstractODEProblem <: DEProblem
  abstract AbstractSDEProblem <: DEProblem
  abstract AbstractDAEProblem <: DEProblem
  abstract AbstractDDEProblem <: DEProblem
  abstract AbstractPoissonProblem <: DEProblem
  abstract AbstractHeatProblem <: DEProblem
  "`PdeSolution`: Wrapper for the objects obtained from a solver"
  abstract DESolution
  abstract AbstractODESolution <: DESolution
  abstract AbstractSDESolution <: AbstractODESolution # Needed for plot recipes
  abstract AbstractFEMSolution <: DESolution
  abstract DESensitivity
  "`Mesh`: An abstract type which holds a (node,elem) pair and other information for a mesh"
  abstract Mesh
  "`Tableau`: Holds the information for a Runge-Kutta Tableau"
  abstract Tableau

  "`ODERKTableau`: A Runge-Kutta Tableau for an ODE integrator"
  abstract ODERKTableau <: Tableau

  "`DEIntegrator`: A DifferentialEquations Integrator type, used to initiate a solver."
  abstract DEIntegrator
  "`DEParameters`: Holds the parameters used in a DifferntialEquations model"
  abstract DEParameters

  include("utils.jl")
  include("solutions.jl")
  include("plotrecipes.jl")
  include("tableaus.jl")

  function solve end

  export DEProblem, DESolution, DEParameters, AbstractDAEProblem, AbstractDDEProblem,
         AbstractODEProblem, AbstractSDEProblem, DAESolution, DEIntegrator, Mesh,
         Tableau, DESensitivity, AbstractODESolution, AbstractPoissonProblem,
         AbstractHeatProblem, AbstractFEMSolution, ODERKTableau, ExplicitRKTableau,
         ImplicitRKTableau, AbstractSDESolution

  export recursivecopy!, @def, vecvecapply, numparameters, copyat_or_push!

end # module
