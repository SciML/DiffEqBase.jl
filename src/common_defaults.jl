function abs2_and_sum(x, y)
    reduce(Base.add_sum, x, init = zero(real(value(eltype(x))))) +
    reduce(Base.add_sum, y, init = zero(real(value(eltype(y)))))
end
@inline UNITLESS_ABS2(x::Number) = abs2(x)
@inline function UNITLESS_ABS2(x::AbstractArray)
    mapreduce(UNITLESS_ABS2, abs2_and_sum, x, init = zero(real(value(eltype(x)))))
end
@inline function UNITLESS_ABS2(x::RecursiveArrayTools.ArrayPartition)
    mapreduce(UNITLESS_ABS2, abs2_and_sum, x.x, init = zero(real(value(eltype(x)))))
end

@inline recursive_length(u::AbstractArray{<:Number}) = length(u)
@inline recursive_length(u::Number) = length(u)
@inline recursive_length(u::AbstractArray{<:AbstractArray}) = sum(recursive_length, u)
@inline recursive_length(u::RecursiveArrayTools.ArrayPartition) = sum(recursive_length, u.x)
@inline recursive_length(u::RecursiveArrayTools.VectorOfArray) = sum(recursive_length, u.u)
@inline function recursive_length(u::AbstractArray{<:StaticArray{S, <:Number}}) where {S}
    prod(Size(eltype(u))) * length(u)
end

@inline ODE_DEFAULT_NORM(u::Union{AbstractFloat, Complex}, t) = @fastmath abs(u)

@inline function ODE_DEFAULT_NORM(u::Array{T}, t) where {T <: Union{AbstractFloat, Complex}}
    x = zero(T)
    @inbounds @fastmath for ui in u
        x += abs2(ui)
    end
    Base.FastMath.sqrt_fast(real(x) / max(length(u), 1))
end

@inline function ODE_DEFAULT_NORM(u::StaticArrays.StaticArray{<:Tuple, T},
                                  t) where {T <: Union{AbstractFloat, Complex}}
    Base.FastMath.sqrt_fast(real(sum(abs2, u)) / max(length(u), 1))
end

@inline function ODE_DEFAULT_NORM(u::AbstractArray, t)
    Base.FastMath.sqrt_fast(UNITLESS_ABS2(u) / max(recursive_length(u), 1))
end
@inline ODE_DEFAULT_NORM(u, t) = norm(u)

@inline ODE_DEFAULT_ISOUTOFDOMAIN(u, p, t) = false
@inline function ODE_DEFAULT_PROG_MESSAGE(dt, u::Array, p, t)
    tmp = u[1]
    for i in eachindex(u)
        tmp = ifelse(abs(u[i]) > abs(tmp), u[i], tmp)
    end
    "dt=" * string(dt) * "\nt=" * string(t) * "\nmax u=" * string(tmp)
end
@inline function ODE_DEFAULT_PROG_MESSAGE(dt, u, p, t)
    "dt=" * string(dt) * "\nt=" * string(t) * "\nmax u=" * string(maximum(abs.(u)))
end

@inline NAN_CHECK(x::Number) = isnan(x)
@inline NAN_CHECK(x::Float64) = isnan(x) || (x > 1e50)
@inline NAN_CHECK(x::Enum) = false
@inline NAN_CHECK(x::AbstractArray) = any(NAN_CHECK, x)
@inline NAN_CHECK(x::RecursiveArrayTools.ArrayPartition) = any(NAN_CHECK, x.x)
@inline function NAN_CHECK(x::SparseArrays.AbstractSparseMatrixCSC)
    any(NAN_CHECK, SparseArrays.nonzeros(x))
end

@inline ODE_DEFAULT_UNSTABLE_CHECK(dt, u, p, t) = false
@inline ODE_DEFAULT_UNSTABLE_CHECK(dt, u::Union{Number, AbstractArray}, p, t) = NAN_CHECK(u)
