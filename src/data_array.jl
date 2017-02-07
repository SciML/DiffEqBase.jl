similar(A::DEDataArray) = deepcopy(A)

done(A::DEDataArray, i::Int64) = done(A.x,i)
eachindex(A::DEDataArray)      = eachindex(A.x)
next(A::DEDataArray, i::Int64) = next(A.x,i)
start(A::DEDataArray)          = start(A.x)

length(A::DEDataArray) = length(A.x)
ndims(A::DEDataArray)  = ndims(A.x)
size(A::DEDataArray)   = size(A.x)

@generated function recursivecopy!{T}(dest::T, src::T)
   assignments = [:(typeof(getfield(dest,$i)) <: AbstractArray ? recursivecopy!(getfield(dest, $i),getfield(src, $i)) : setfield!(dest, $i, getfield(src, $i))) for i=1:nfields(T)]
   :($(assignments...); dest)
end

getindex( A::DEDataArray,    i::Int) = (A.x[i])
setindex!(A::DEDataArray, x, i::Int) = (A.x[i] = x)
getindex( A::DEDataArray,    i::Int...) = (A.x[i...])
setindex!(A::DEDataArray, x, i::Int...) = (A.x[i...] = x)
Base.linearindexing{T<:DEDataArray}(::Type{T}) = Base.LinearFast()
