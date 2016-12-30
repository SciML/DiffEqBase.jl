"""
`ExplicitRKTableau`

Holds a tableau which defines an explicit Runge-Kutta method.
"""
type ExplicitRKTableau{MType<:AbstractMatrix,VType<:AbstractVector} <: ODERKTableau
  A::MType
  c::VType
  α::VType
  αEEst::VType
  stages::Int
  order::Int
  adaptiveorder::Int #The lower order of the pair. Only used for adaptivity.
  fsal::Bool # First same as last
end
ExplicitRKTableau(A,c,α,order;adaptiveorder=0,αEEst=Float64[],fsal=false) = ExplicitRKTableau(A,c,α,αEEst,length(α),order,adaptiveorder,fsal)

"""
`ImplicitRKTableau`

Holds a tableau which defines an implicit Runge-Kutta method.
"""
type ImplicitRKTableau{MType<:AbstractMatrix,VType<:AbstractVector} <: ODERKTableau
  A::MType
  c::VType
  α::VType
  αEEst::VType
  stages::Int
  order::Int
  adaptiveorder::Int #The lower order of the pair. Only used for adaptivity.
  ImplicitRKTableau(A,c,α,order;adaptiveorder=0,αEEst=Float64[]) = new(A,c,α,αEEst,length(α),order,adaptiveorder)
end
ImplicitRKTableau(A,c,α,order;adaptiveorder=0,αEEst=Float64[],fsal=false) = ImplicitRKTableau(A,c,α,αEEst,length(α),order,adaptiveorder,fsal)
