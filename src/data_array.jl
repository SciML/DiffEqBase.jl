similar(A::DEDataArray) = deepcopy(A)

@generated function similar(A::DEDataArray,dims::Tuple)
  assignments = [s == :x ? :(similar(A.x,dims)) : (s_new = Meta.quot(:($s)); :(deepcopy(getfield(A,$s_new)))) for s in fieldnames(A)]
  :(typeof(A)($(assignments...)))
end

@generated function similar{T,N}(A::DEDataArray,::Type{T},dims::Tuple{Vararg{Int,N}})
  assignments = [s == :x ? :(zeros(A.x,T,dims)) : (s_new = Meta.quot(:($s)); :(deepcopy(getfield(A,$s_new)))) for s in fieldnames(A)]
  :(parameterless_type(A)($(assignments...)))
end

done(A::DEDataArray, i::Int) = done(A.x,i)
eachindex(A::DEDataArray)      = eachindex(A.x)
next(A::DEDataArray, i::Int) = next(A.x,i)
start(A::DEDataArray)          = start(A.x)

length(A::DEDataArray) = length(A.x)
ndims(A::DEDataArray)  = ndims(A.x)
size(A::DEDataArray)   = size(A.x)

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

getindex( A::DEDataArray,    i::Int) = (A.x[i])
setindex!(A::DEDataArray, x, i::Int) = (A.x[i] = x)
getindex( A::DEDataArray,    i::Int...) = (A.x[i...])
setindex!(A::DEDataArray, x, i::Int...) = (A.x[i...] = x)
Base.IndexStyle(::Type{DEDataArray}) = Base.IndexLinear()
