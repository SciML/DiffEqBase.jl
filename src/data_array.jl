similar(A::DEDataArray) = deepcopy(A)

done(A::DEDataArray, i::Int64) = done(A.x,i)
eachindex(A::DEDataArray)      = eachindex(A.x)
next(A::DEDataArray, i::Int64) = next(A.x,i)
start(A::DEDataArray)          = start(A.x)

length(A::DEDataArray) = length(A.x)
ndims(A::DEDataArray)  = ndims(A.x)
size(A::DEDataArray)   = size(A.x)

# Maybe should use fieldtype(typeof(dest), $i) ?
@generated function recursivecopy!{T<:DEDataArray}(dest::T, src::T)
   assignments = [:(typeof(getfield(dest,$i)) <: AbstractArray ? recursivecopy!(getfield(dest, $i),getfield(src, $i)) : setfield!(dest, $i, getfield(src, $i))) for i=1:nfields(T)]
   :($(assignments...); dest)
end

@generated function copy_non_array_fields{T<:DEDataArray}(previous::T, arr::AbstractArray)
  assignments = [:(getfield(previous,$i)) for i=2:nfields(T)]
  :(T(arr,$(assignments...)))
end

getindex( A::DEDataArray,    i::Int) = (A.x[i])
setindex!(A::DEDataArray, x, i::Int) = (A.x[i] = x)
getindex( A::DEDataArray,    i::Int...) = (A.x[i...])
setindex!(A::DEDataArray, x, i::Int...) = (A.x[i...] = x)
@compat Base.IndexStyle(::Type{DEDataArray}) = Base.IndexLinear()
