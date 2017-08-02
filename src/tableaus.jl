"""
`ExplicitRKTableau`

Holds a tableau which defines an explicit Runge-Kutta method.
"""
mutable struct ExplicitRKTableau{MType<:AbstractMatrix,VType<:AbstractVector,fsal} <: ODERKTableau
  A::MType
  c::VType
  α::VType
  αEEst::VType
  stages::Int
  order::Int
  adaptiveorder::Int #The lower order of the pair. Only used for adaptivity.
end
ExplicitRKTableau(A::MType,c::VType,α::VType,order;
                  adaptiveorder=0,αEEst=VType(),fsal=false) where {MType,VType} =
                  ExplicitRKTableau{MType,VType,fsal}(
                  A,c,α,αEEst,length(α),order,adaptiveorder)

"""
`ImplicitRKTableau`

Holds a tableau which defines an implicit Runge-Kutta method.
"""
mutable struct ImplicitRKTableau{MType<:AbstractMatrix,VType<:AbstractVector} <: ODERKTableau
  A::MType
  c::VType
  α::VType
  αEEst::VType
  stages::Int
  order::Int
  adaptiveorder::Int #The lower order of the pair. Only used for adaptivity.
end
ImplicitRKTableau(A::MType,c::VType,α::VType,order;
                  adaptiveorder=0,αEEst=VType()) where {MType,VType} =
                  ImplicitRKTableau{MType,VType}(
                  A,c,α,αEEst,length(α),order,adaptiveorder)
