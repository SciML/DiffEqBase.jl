function recursivecopy!{T<:Number,N}(b::Array{T,N},a::Array{T,N})
  @inbounds copy!(b,a)
end

function recursivecopy!{T<:AbstractArray,N}(b::Array{T,N},a::Array{T,N})
  @inbounds for i in eachindex(a)
    recursivecopy!(b[i],a[i])
  end
end

macro def(name, definition)
    quote
        macro $name()
            esc($(Expr(:quote, definition)))
        end
    end
end


function vecvecapply{T<:Number,N}(f::Base.Callable,v::Vector{Array{T,N}})
  sol = Vector{eltype(eltype(v))}(0)
  for i in eachindex(v)
    for j in eachindex(v[i])
      push!(sol,v[i][j])
    end
  end
  f(sol)
end

function vecvecapply{T<:Number}(f::Base.Callable,v::Vector{T})
  f(v)
end

"""
`numparameters(f)`

Returns the number of parameters of `f` for the method which has the most parameters.
"""
function numparameters(f)
  #=
  if length(methods(f))>1
    warn("Number of methods for f is greater than 1. Choosing linearity based off of method with most parameters")
  end
  =#
  numparm = maximum([length(m.sig.parameters) for m in methods(f)]) #in v0.5, all are generic
  return (numparm-1) #-1 in v0.5 since it adds f as the first parameter
end

@inline function copyat_or_push!{T}(a::AbstractVector{T},i::Int,x)
  @inbounds if length(a) >= i
    if T <: Number
      a[i] = x
    else
      recursivecopy!(a[i],x)
    end
  else
    if T <: Number
      push!(a,copy(x))
    else
      push!(a,deepcopy(x))
    end
  end
  nothing
end
