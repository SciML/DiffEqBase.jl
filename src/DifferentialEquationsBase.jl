module DifferentialEquationsBase

  "`DEProblem`: Defines differential equation problems via its internal functions"
  abstract DEProblem
  abstract DEElement
  abstract AbstractODEProblem <: DEProblem
  abstract AbstractSDEProblem <: DEProblem
  abstract AbstractDAEProblem <: DEProblem
  abstract AbstractDDEProblem <: DEProblem
  "`PdeSolution`: Wrapper for the objects obtained from a solver"
  abstract DESolution
  abstract AbstractODESolution <: DESolution
  abstract DESensitivity
  "`Mesh`: An abstract type which holds a (node,elem) pair and other information for a mesh"
  abstract Mesh
  "`Tableau`: Holds the information for a Runge-Kutta Tableau"
  abstract Tableau
  "`DEIntegrator`: A DifferentialEquations Integrator type, used to initiate a solver."
  abstract DEIntegrator
  "`DEParameters`: Holds the parameters used in a DifferntialEquations model"
  abstract DEParameters




  function recursivecopy!{T<:Number,N}(b::Array{T,N},a::Array{T,N})
    @inbounds copy!(b,a)
  end

  function recursivecopy!{T<:AbstractArray,N}(b::Array{T,N},a::Array{T,N})
    @inbounds for i in eachindex(a)
      recursivecopy!(b[i],a[i])
    end
  end

  export DEProblem, DESolution, DEParameters, AbstractDAEProblem, AbstractDDEProblem,
         AbstractODEProblem, AbstractSDEProblem, DAESolution, DEIntegrator, Mesh,
         Tableau, DESensitivity, AbstractODESolution

  export recursivecopy!


end # module
