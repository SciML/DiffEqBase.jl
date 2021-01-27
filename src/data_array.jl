# iteration
Base.iterate(A::DEDataArray) = iterate(A.x)
Base.iterate(A::DEDataArray, state) = iterate(A.x, state)

# size
Base.length(A::DEDataArray) = length(A.x)
Base.size(A::DEDataArray) = size(A.x)

# indexing
@inline function Base.getindex(A::DEDataArray, I...)
    @boundscheck checkbounds(A.x, I...)
    @inbounds return A.x[I...]
end
@inline function Base.setindex!(A::DEDataArray, x, I...)
    @boundscheck checkbounds(A.x, I...)
    @inbounds A.x[I...] = x
end
Base.axes(A::DEDataArray) = axes(A.x)
Base.LinearIndices(A::DEDataArray) = LinearIndices(A.x)
Base.IndexStyle(::Type{<:DEDataArray}) = Base.IndexLinear()

Base.copy(A::DEDataArray) = deepcopy(A)

# zero data arrays
@generated function Base.zero(A::DEDataArray)
    assignments = [s == :x ? :(zero(A.x)) :
                   (sq = Meta.quot(s); :(deepcopy(getfield(A, $sq))))
                   for s in fieldnames(A)]
    :(DiffEqBase.parameterless_type(A)($(assignments...)))
end

# similar data arrays
@generated function Base.similar(A::DEDataArray, ::Type{T}, dims::NTuple{N,Int}) where {T,N}
    assignments = [s == :x ? :(typeof(A.x) <: StaticArray ? similar(A.x, T, Size(A.x)) : similar(A.x, T, dims)) :
                   (sq = Meta.quot(s); :(deepcopy(getfield(A, $sq))))
                   for s in fieldnames(A)]
    :(DiffEqBase.parameterless_type(A)($(assignments...)))
end

@generated function Base.similar(A::DEDataArray, ::Type{T}) where {T}
    assignments = [s == :x ? :(typeof(A.x) <: StaticArray ? similar(A.x, T, Size(A.x)) : similar(A.x, T)) :
                   (sq = Meta.quot(s); :(deepcopy(getfield(A, $sq))))
                   for s in fieldnames(A)]
    :(DiffEqBase.parameterless_type(A)($(assignments...)))
end

@generated function Base.similar(A::DEDataArray) where {T}
    assignments = [s == :x ? :(typeof(A.x) <: StaticArray ? similar(A.x) : similar(A.x)) :
                   (sq = Meta.quot(s); :(deepcopy(getfield(A, $sq))))
                   for s in fieldnames(A)]
    :(DiffEqBase.parameterless_type(A)($(assignments...)))
end

"""
    recursivecopy!(dest::T, src::T) where {T<:DEDataArray}

Recursively copy fields of `src` to `dest`.
"""
@generated function RecursiveArrayTools.recursivecopy!(dest::T, src::T) where {T<:DEDataArray}
    fields = fieldnames(src)

    expressions = Vector{Expr}(undef, length(fields))

    @inbounds for i = 1:length(fields)
        f  = fields[i]
        Tf = src.types[i]
        qf = Meta.quot(f)

        if !ArrayInterface.ismutable(Tf)
            expressions[i] = :( dest.$f = getfield( src, $qf ) )
        elseif Tf <: AbstractArray
            expressions[i] = :( RecursiveArrayTools.recursivecopy!(dest.$f, getfield( src, $qf ) ) )
        else
            expressions[i] = :( dest.$f = deepcopy( getfield( src, $qf ) ) )
        end
    end

    :($(expressions...); dest)
end

"""
    copy_fields(arr:AbstractArray, template::DEDataArray)

Create `DEDataArray` that wraps `arr` with all other fields set to a deep copy of the
value in `template`.

copy_fields!(dest::T, src::T2) where {T<:DEDataArray,T2<:DEDataArray}

Replace all fields of `dest` except of its wrapped array with a copy of the
value in `src`. Arrays are recursively copied.
"""
@generated function copy_fields(arr::AbstractArray, template::DEDataArray)
    assignments = [s == :x ? :(arr) :
                   (sq = Meta.quot(s); :(deepcopy(getfield(template, $sq))))
                   for s in fieldnames(template)]
    :(parameterless_type(template)($(assignments...)))
end

@generated function copy_fields!(dest::T, src::T2) where
    {T<:DEDataArray,T2<:DEDataArray}

    fields = fieldnames(src)

    expressions = Vector{Expr}(undef, length(fields))

    @inbounds for i = 1:length(fields)
        f  = fields[i]
        Tf = src.types[i]
        qf = Meta.quot(f)

        if f == :x
            expressions[i] = :( )
        elseif !ArrayInterface.ismutable(Tf)
            expressions[i] = :( dest.$f = getfield( src, $qf ) )
        elseif Tf <: AbstractArray
            expressions[i] = :( RecursiveArrayTools.recursivecopy!(dest.$f, getfield( src, $qf ) ) )
        else
            expressions[i] = :( dest.$f = deepcopy( getfield( src, $qf ) ) )
        end
    end

    :($(expressions...); dest)
end

################# Overloads for stiff solvers ##################################

LinearAlgebra.ldiv!(A::DEDataArray,F::Factorization, B::DEDataArray) = ldiv!(A.x,F,B.x)
LinearAlgebra.ldiv!(F::Factorization, B::DEDataArray) = ldiv!(F, B.x)
LinearAlgebra.ldiv!(F::Factorization,A::Base.ReshapedArray{T1,T2,T3,T4}) where {T1,T2,T3<:DEDataArray,T4} = ldiv!(F,vec(A.parent.x))
Base.:+(::LinearAlgebra.UniformScaling,x::DEDataArray) = DiffEqBase.copy_fields(I + x.x,x)

Base.unsafe_convert(::Type{Ptr{T}}, a::DEDataArray{T}) where {T} = Base.unsafe_convert(Ptr{T}, getfield(a,:x))
ArrayInterface.zeromatrix(x::DEDataArray) = ArrayInterface.zeromatrix(x.x)

################# Broadcast ####################################################

const DEDataArrayStyle = Broadcast.ArrayStyle{DEDataArray}
Base.BroadcastStyle(::Type{<:DEDataArray}) = Broadcast.ArrayStyle{DEDataArray}()
Base.BroadcastStyle(::Broadcast.ArrayStyle{DEDataArray},::Broadcast.ArrayStyle) = Broadcast.ArrayStyle{DEDataArray}()
Base.BroadcastStyle(::Broadcast.ArrayStyle,::Broadcast.ArrayStyle{DEDataArray}) = Broadcast.ArrayStyle{DEDataArray}()
Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{DEDataArray}},::Type{ElType}) where ElType = similar(find_dedata(bc))

find_dedata(bc::Base.Broadcast.Broadcasted) = find_dedata(bc.args)
function find_dedata(args::Tuple)
  !isempty(args) && find_dedata(find_dedata(args[1]), Base.tail(args))
end
find_dedata(x) = x
find_dedata(a::DEDataArray, rest) = a
find_dedata(::Any, rest) = find_dedata(rest)

@inline function Base.copy(bc::Broadcast.Broadcasted{DEDataArrayStyle})
    out = find_dedata(bc)
    copy_fields(copy(unpack(bc)), out)
end

# drop DEData part
@inline unpack(bc::Broadcast.Broadcasted{Style}) where Style = Broadcast.Broadcasted{Style}(bc.f, unpack_args(bc.args))
@inline unpack(bc::Broadcast.Broadcasted{DEDataArrayStyle}) = Broadcast.Broadcasted(bc.f, unpack_args(bc.args))
unpack(x) = x
unpack(x::DEDataArray) = x.x

@inline unpack_args(args::Tuple) = (unpack(args[1]), unpack_args(Base.tail(args))...)
unpack_args(args::Tuple{Any}) = (unpack(args[1]),)
unpack_args(::Any, args::Tuple{}) = ()

# Broadcasting checks for aliasing with Base.dataids but the fallback
# for AbstractArrays is very slow. Instead, we just call dataids on the
# wrapped buffer
Base.dataids(A::DEDataArray) = Base.dataids(A.x)
