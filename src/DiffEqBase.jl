__precompile__()

module DiffEqBase

  using RecipesBase, Parameters, RecursiveArrayTools, SimpleTraits
  using Ranges # For plot recipes with units
  import Base: length, size, getindex, endof, show, print,
               next, start, done, eltype, eachindex, convert, promote

  ## Problems
  "`DEProblem`: Defines differential equation problems via its internal functions"
  abstract DEProblem
  abstract DEElement
  abstract DESensitivity

  # Initial Value Problems (IVP)
  abstract AbstractIVPProblem{uType,tType,isinplace,F} <: DEProblem
  "Checks whether uses an in-place function or not"
  isinplace{uType,tType,I,F}(::Type{AbstractIVPProblem{uType,tType,I,F}}) = I
  isinplace{A<:AbstractIVPProblem}(::Type{A}) = isinplace(supertype(A))
  isinplace{A<:AbstractIVPProblem}(::A) = isinplace(A)
  @traitdef IsInplace{P}
  @traitimpl IsInplace{P} <- isinplace(P)

  # IVP subtypes
  abstract AbstractODEProblem{uType,tType,isinplace,F} <: AbstractIVPProblem{uType,tType,isinplace,F}
  abstract AbstractODETestProblem{uType,tType,isinplace,F} <: AbstractODEProblem{uType,tType,isinplace,F}
  abstract AbstractSDEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3}  <: AbstractIVPProblem{uType,tType,isinplace,F}
  abstract AbstractSDETestProblem{uType,tType,isinplace,NoiseClass,F,F2,F3}  <: AbstractSDEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3}
  abstract AbstractDAEProblem{uType,duType,tType,isinplace,F} <: AbstractIVPProblem{uType,tType,isinplace,F}
  abstract AbstractDAETestProblem{uType,duType,tType,isinplace,F} <: AbstractDAEProblem{uType,duType,tType,isinplace,F}

  # Misc other problems
  abstract AbstractDDEProblem <: DEProblem
  abstract AbstractPoissonProblem <: DEProblem
  abstract AbstractHeatProblem <: DEProblem

  ## Algorithms
  abstract DEAlgorithm
  abstract IVPAlgorithm <: DEAlgorithm
  abstract AbstractODEAlgorithm <: IVPAlgorithm
  abstract AbstractSDEAlgorithm <: IVPAlgorithm
  abstract AbstractDAEAlgorithm <: IVPAlgorithm

  ## Solutions
  abstract DESolution
  abstract AbstractODESolution <: DESolution
  abstract AbstractSDESolution <: AbstractODESolution # Needed for plot recipes
  abstract AbstractDAESolution <: AbstractODESolution # Needed for plot recipes
  abstract AbstractDDESolution <: AbstractODESolution # Needed for plot recipes
  abstract AbstractFEMSolution <: DESolution
  abstract AbstractSensitivitySolution

  ## Misc
  "`Mesh`: An abstract type which holds a (node,elem) pair and other information for a mesh"
  abstract Mesh
  "`Tableau`: Holds the information for a Runge-Kutta Tableau"
  abstract Tableau
  "`ODERKTableau`: A Runge-Kutta Tableau for an ODE integrator"
  abstract ODERKTableau <: Tableau

  include("utils.jl")
  include("extended_functions.jl")
  include("noise_process.jl")
  include("solutions/ode_solutions.jl")
  include("solutions/sde_solutions.jl")
  include("solutions/dae_solutions.jl")
  include("solutions/solution_interface.jl")
  include("tableaus.jl")
  include("problems/ode_problems.jl")
  include("problems/sde_problems.jl")
  include("problems/dae_problems.jl")

  solve(p::DEProblem, a::DEAlgorithm; kwargs...) = solve(promote(p, a), a; kwargs...)
  promote(p::DEProblem, a::DEAlgorithm) =
      error("Algorithm (solver) of type `$(typeof(a))` is not applicable to problem of type `$(typeof(p))`")
  convert{P1<:DEProblem}(::Type{P1}, p::DEProblem) =
      error("Cannot convert problem type $(typeof(p)) into problem type $P1")
  convert{P1<:DEProblem}(::Type{P1}, p::P1) = p

  export DEProblem, DESolution, DEParameters, AbstractDAEProblem, AbstractDDEProblem,
         AbstractODEProblem, AbstractSDEProblem, DAESolution, DEIntegrator, Mesh,
         Tableau, DESensitivity, AbstractODESolution, AbstractPoissonProblem,
         AbstractHeatProblem, AbstractFEMSolution, ODERKTableau, ExplicitRKTableau,
         ImplicitRKTableau, AbstractSDESolution, DESensitivity, solve

  export tuples

  export numparameters

  export HasJac, HastGrad, HasParamFuncs, HasParamDeriv, HasParamJac, HasInvJac,
         HasInvW, HasInvW_t, HasHes, HasInvHes

  export has_jac, has_invjac, has_invW, has_invW_t, has_hes, has_invhes,
         has_tgrad, has_paramfuncs, has_paramderiv, has_paramjac

  export ODEProblem, ODETestProblem, ODESolution, ODETestSolution

  export SDEProblem, SDETestProblem, SDESolution, SDETestSolution

  export DAEProblem, DAETestProblem, DAESolution, DAETestSolution

  export build_solution

  # Algorithms

  export AbstractODEAlgorithm, AbstractSDEAlgorithm, AbstractDAEAlgorithm

end # module
