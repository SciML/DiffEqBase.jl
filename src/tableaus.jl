# Maybe remove these?

"""
`ExplicitRKTableau`

Holds a tableau which defines an explicit Runge-Kutta method.
"""
type ExplicitRKTableau <: ODERKTableau
  A#::Array{Float64,2}
  c#::Vector{Float64}
  α#::Vector{Float64}
  αEEst#::Vector{Float64}
  stages::Int
  order::Int
  adaptiveorder::Int #The lower order of the pair. Only used for adaptivity.
  fsal::Bool # First same as last
  ExplicitRKTableau(A,c,α,order;adaptiveorder=0,αEEst=Float64[],fsal=false) = new(A,c,α,αEEst,length(α),order,adaptiveorder,fsal)
end

"""
`ImplicitRKTableau`

Holds a tableau which defines an implicit Runge-Kutta method.
"""
type ImplicitRKTableau <: ODERKTableau
  A#::Array{Float64,2}
  c#::Vector{Float64}
  α#::Vector{Float64}
  αEEst#::Vector{Float64}
  stages::Int
  order::Int
  adaptiveorder::Int #The lower order of the pair. Only used for adaptivity.
  ImplicitRKTableau(A,c,α,order;adaptiveorder=0,αEEst=Float64[]) = new(A,c,α,αEEst,length(α),order,adaptiveorder)
end
