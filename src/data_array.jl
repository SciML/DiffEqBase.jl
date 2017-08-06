# define AbstractArray interface of DEDataArray

# iteration
start(A::DEDataArray)  = start(A.x)
next(A::DEDataArray, i) = next(A.x, i)
done(A::DEDataArray, i) = done(A.x, i)
length(A::DEDataArray) = length(A.x)

# size
size(A::DEDataArray) = (length(A.x),)
ndims(A::DEDataArray) = 1

# indexing
getindex(A::DEDataArray, i::Int) = (A.x[i])
setindex!(A::DEDataArray, x, i::Int) = (A.x[i] = x)
Base.IndexStyle(::Type{<:DEDataArray}) = Base.IndexLinear()

similar(A::DEDataArray) = deepcopy(A)

@generated function similar(A::DEDataArray,dims::Tuple)
  assignments = [s == :x ? :(similar(A.x,dims)) : (s_new = Meta.quot(:($s)); :(deepcopy(getfield(A,$s_new)))) for s in fieldnames(A)]
  :(typeof(A)($(assignments...)))
end

@generated function similar{T,N}(A::DEDataArray,::Type{T},dims::Tuple{Vararg{Int,N}})
  assignments = [s == :x ? :(zeros(A.x,T,dims)) : (s_new = Meta.quot(:($s)); :(deepcopy(getfield(A,$s_new)))) for s in fieldnames(A)]
  :(parameterless_type(A)($(assignments...)))
end

# Maybe should use fieldtype(typeof(dest), $i) ?
@generated function recursivecopy!{T}(dest::DEDataArray{T}, src::DEDataArray{T})
   assignments = [:(typeof(getfield(dest,$i)) <: AbstractArray ? recursivecopy!(getfield(dest, $i),getfield(src, $i)) : setfield!(dest, $i, getfield(src, $i))) for i=1:nfields(dest)]
   :($(assignments...); dest)
end

@generated function copy_non_array_fields{T}(arr::AbstractArray, previous::DEDataArray{T})
  assignments = [:(getfield(previous,$i)) for i=2:nfields(previous)]
  :(typeof(previous)(arr,$(assignments...)))
end

@generated function copy_non_array_fields!{T}(dest::DEDataArray{T}, src::DEDataArray{T})
  assignments = [:(typeof(getfield(dest,$i)) <: AbstractArray ? recursivecopy!(getfield(dest, $i),getfield(src, $i)) : setfield!(dest, $i, getfield(src, $i))) for i=2:nfields(dest)]
  :($(assignments...); dest)
end

# ensure that broadcasts dispatch to our broadcast methods below
_containertype(::Type{T}) where {T<:DEDataArray} = T
promote_containertype(::Type{T}, ::Type) where {T<:DEDataArray} = T
promote_containertype(::Type, ::Type{T}) where {T<:DEDataArray} = T
# avoid ambiguous definitions
promote_containertype(::Type{T}, ::Type{T}) where {T<:DEDataArray} = T
promote_containertype(::Type{T}, ::Type{Array}) where {T<:DEDataArray} = T
promote_containertype(::Type{Array}, ::Type{T}) where {T<:DEDataArray} = T

# apply generic broadcast methods to contained arrays
broadcast_c(f, ::Type{<:DEDataArray}, as...) = broadcast(f, map(broadcast_content, as)...)

# update contained array with generic broadcast methods
function broadcast_c!(f, ::Type{<:DEDataArray}, ::Type, dest, as...)
    broadcast!(f, dest.x, map(broadcast_content, as)...)
    dest
end

broadcast_content(x) = x
broadcast_content(x::DEDataArray) = x.x
