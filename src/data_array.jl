# iteration
start(A::DEDataArray) = start(A.x)
next(A::DEDataArray, i) = next(A.x, i)
done(A::DEDataArray, i) = done(A.x, i)

# size
length(A::DEDataArray) = length(A.x)
size(A::DEDataArray) = size(A.x)

# indexing
@inline function getindex(A::DEDataArray, I...)
    @boundscheck checkbounds(A.x, I...)
    @inbounds return A.x[I...]
end
@inline function setindex!(A::DEDataArray, x, I...)
    @boundscheck checkbounds(A.x, I...)
    @inbounds A.x[I...] = x
end
indices(A::DEDataArray) = indices(A.x)
Base.LinearIndices(A::DEDataArray) = LinearIndices(A.x)
Base.IndexStyle(::Type{<:DEDataArray}) = Base.IndexLinear()

# similar data arrays
@generated function similar(A::DEDataArray, ::Type{T}, dims::NTuple{N,Int}) where {T,N}
    assignments = [s == :x ? :(similar(A.x, T, dims)) :
                   (sq = Meta.quot(s); :(deepcopy(getfield(A, $sq))))
                   for s in fieldnames(A)]
    :(parameterless_type(A)($(assignments...)))
end

"""
    recursivecopy!(dest::T, src::T) where {T<:DEDataArray}

Recursively copy fields of `src` to `dest`.
"""
@generated function recursivecopy!(dest::T, src::T) where {T<:DEDataArray}
    assignments = [(sq = Meta.quot(s);
                   :(typeof(getfield(dest, $sq)) <: AbstractArray ?
                     recursivecopy!(getfield(dest, $sq), getfield(src, $sq)) :
                     setfield!(dest, $sq, deepcopy(getfield(src, $sq)))))
                   for s in fieldnames(dest)]
    :($(assignments...); dest)
end

"""
    copy_fields(arr:AbstractArray, template::DEDataArray)

Create `DEDataArray` that wraps `arr` with all other fields set to a deep copy of the
value in `template`.
"""
@generated function copy_fields(arr::AbstractArray, template::DEDataArray)
    assignments = [s == :x ? :(arr) :
                   (sq = Meta.quot(s); :(deepcopy(getfield(template, $sq))))
                   for s in fieldnames(template)]
    :(parameterless_type(template)($(assignments...)))
end

"""
    copy_fields!(dest::T, src::T) where {T<:DEDataArray}

Replace all fields of `dest` except of its wrapped array with a copy of the value in `src`.
Arrays are recursively copied.
"""
@generated function copy_fields!(dest::T, src::T2) where
    {T<:DEDataArray,T2<:DEDataArray}
    assignments = [(sq = Meta.quot(s);
                    :(typeof(getfield(dest, $sq)) <: AbstractArray ?
                      recursivecopy!(getfield(dest, $sq), getfield(src, $sq)) :
                      setfield!(dest, $sq, deepcopy(getfield(src, $sq)))))
                   for s in fieldnames(dest) if s != :x]
    :($(assignments...); dest)
end

################# Overloads for stiff solvers ##################################

LinearAlgebra.ldiv!(A::DEDataArray,F::Factorization, B::DEDataArray) = ldiv!(A.x,F,B.x)
LinearAlgebra.ldiv!(F::Factorization, B::DEDataArray) = ldiv!(F, B.x)
LinearAlgebra.ldiv!(F::Factorization,A::Base.ReshapedArray{T1,T2,T3,T4}) where {T1,T2,T3<:DEDataArray,T4} = ldiv!(F,vec(A.parent.x))
