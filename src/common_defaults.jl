function abs2_and_sum(x, y)
    reduce(Base.add_sum, x, init = zero(real(value(eltype(x))))) +
    reduce(Base.add_sum, y, init = zero(real(value(eltype(y)))))
end
@inline UNITLESS_ABS2(x::Number) = abs2(x)
@inline function UNITLESS_ABS2(x::AbstractArray)
    mapreduce(UNITLESS_ABS2, abs2_and_sum, x, init = zero(real(value(eltype(x)))))
end
@inline function UNITLESS_ABS2(x::RecursiveArrayTools.AbstractVectorOfArray)
    mapreduce(UNITLESS_ABS2, abs2_and_sum, x.u, init = zero(real(value(eltype(x)))))
end
@inline function UNITLESS_ABS2(x::RecursiveArrayTools.ArrayPartition)
    mapreduce(UNITLESS_ABS2, abs2_and_sum, x.x, init = zero(real(value(eltype(x)))))
end

@inline UNITLESS_ABS2(f::F, x::Number) where {F} = abs2(f(x))
@inline function UNITLESS_ABS2(f::F, x::AbstractArray) where {F}
    return mapreduce(UNITLESS_ABS2 ∘ f, abs2_and_sum, x;
        init = zero(real(value(eltype(x)))))
end
@inline function UNITLESS_ABS2(f::F, x::RecursiveArrayTools.ArrayPartition) where {F}
    return mapreduce(UNITLESS_ABS2 ∘ f, abs2_and_sum, x.x;
        init = zero(real(value(eltype(x)))))
end

@inline recursive_length(u::AbstractArray{<:Number}) = length(u)
@inline recursive_length(u::Number) = length(u)
@inline recursive_length(u::AbstractArray{<:AbstractArray}) = sum(recursive_length, u)
@inline recursive_length(u::RecursiveArrayTools.ArrayPartition) = sum(recursive_length, u.x)
@inline recursive_length(u::RecursiveArrayTools.VectorOfArray) = sum(recursive_length, u.u)
@inline function recursive_length(u::AbstractArray{
        <:StaticArraysCore.StaticArray{S, <:Number}}) where {S}
    prod(Size(eltype(u))) * length(u)
end

@inline ODE_DEFAULT_NORM(u::Union{AbstractFloat, Complex}, t) = @fastmath abs(u)

@inline function ODE_DEFAULT_NORM(f::F, u::Union{AbstractFloat, Complex}, t) where {F}
    return @fastmath abs(f(u))
end

@inline function ODE_DEFAULT_NORM(u::Array{T}, t) where {T <: Union{AbstractFloat, Complex}}
    x = zero(T)
    @inbounds @fastmath for ui in u
        x += abs2(ui)
    end
    Base.FastMath.sqrt_fast(real(x) / max(length(u), 1))
end

@inline function ODE_DEFAULT_NORM(f::F,
        u::Union{Array{T}, Iterators.Zip{<:Tuple{Vararg{Array{T}}}}},
        t) where {F, T <: Union{AbstractFloat, Complex}}
    x = zero(T)
    @inbounds @fastmath for ui in u
        x += abs2(f(ui))
    end
    Base.FastMath.sqrt_fast(real(x) / max(length(u), 1))
end

@inline function ODE_DEFAULT_NORM(u::StaticArraysCore.StaticArray{<:Tuple, T},
        t) where {T <: Union{AbstractFloat, Complex}}
    Base.FastMath.sqrt_fast(real(sum(abs2, u)) / max(length(u), 1))
end

@inline function ODE_DEFAULT_NORM(f::F, u::StaticArraysCore.StaticArray{<:Tuple, T},
        t) where {F, T <: Union{AbstractFloat, Complex}}
    Base.FastMath.sqrt_fast(real(sum(abs2 ∘ f, u)) / max(length(u), 1))
end

@inline function ODE_DEFAULT_NORM(
        u::Union{
            AbstractArray,
            RecursiveArrayTools.AbstractVectorOfArray
        },
        t)
    Base.FastMath.sqrt_fast(UNITLESS_ABS2(u) / max(recursive_length(u), 1))
end

@inline function ODE_DEFAULT_NORM(f::F, u::AbstractArray, t) where {F}
    Base.FastMath.sqrt_fast(UNITLESS_ABS2(f, u) / max(recursive_length(u), 1))
end

@inline ODE_DEFAULT_NORM(u, t) = norm(u)
@inline ODE_DEFAULT_NORM(f::F, u, t) where {F} = norm(f.(u))

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

NAN_CHECK(x::Number) = isnan(x)
NAN_CHECK(x::Enum) = false
NAN_CHECK(x::Union{AbstractArray, RecursiveArrayTools.AbstractVectorOfArray}) = any(
    NAN_CHECK, x)
NAN_CHECK(x::RecursiveArrayTools.ArrayPartition) = any(NAN_CHECK, x.x)
function NAN_CHECK(x::SparseArrays.AbstractSparseMatrixCSC)
    any(NAN_CHECK, SparseArrays.nonzeros(x))
end

INFINITE_OR_GIANT(x::Number) = !isfinite(x)
INFINITE_OR_GIANT(x::AbstractFloat) = !(abs(x) < sqrt(floatmax(x)))
INFINITE_OR_GIANT(x::Union{AbstractArray, RecursiveArrayTools.AbstractVectorOfArray}) = any(
    INFINITE_OR_GIANT, x)
INFINITE_OR_GIANT(x::RecursiveArrayTools.ArrayPartition) = any(INFINITE_OR_GIANT, x.x)
function INFINITE_OR_GIANT(x::SparseArrays.AbstractSparseMatrixCSC)
    any(INFINITE_OR_GIANT, SparseArrays.nonzeros(x))
end
ODE_DEFAULT_UNSTABLE_CHECK(dt, u, p, t) = false
function ODE_DEFAULT_UNSTABLE_CHECK(dt, u::Union{Number, AbstractArray{<:Number}}, p, t)
    INFINITE_OR_GIANT(u)
end


# Nonlinear Solve Norm (norm(_, 2))
@inline NONLINEARSOLVE_DEFAULT_NORM(u::Union{AbstractFloat, Complex}) = @fastmath abs(u)
@inline function NONLINEARSOLVE_DEFAULT_NORM(f::F,
        u::Union{AbstractFloat, Complex}) where {F}
    return @fastmath abs(f(u))
end

@inline function NONLINEARSOLVE_DEFAULT_NORM(u::Array{
        T}) where {T <: Union{AbstractFloat, Complex}}
    x = zero(T)
    @inbounds @fastmath for ui in u
        x += abs2(ui)
    end
    return Base.FastMath.sqrt_fast(real(x))
end

@inline function NONLINEARSOLVE_DEFAULT_NORM(f::F,
        u::Union{Array{T}, Iterators.Zip{<:Tuple{Vararg{Array{T}}}}}) where {
        F, T <: Union{AbstractFloat, Complex}}
    x = zero(T)
    @inbounds @fastmath for ui in u
        x += abs2(f(ui))
    end
    return Base.FastMath.sqrt_fast(real(x))
end

@inline function NONLINEARSOLVE_DEFAULT_NORM(u::StaticArraysCore.StaticArray{
        <:Tuple, T}) where {T <: Union{AbstractFloat, Complex}}
    return Base.FastMath.sqrt_fast(real(sum(abs2, u)))
end

@inline function NONLINEARSOLVE_DEFAULT_NORM(f::F,
        u::StaticArraysCore.StaticArray{<:Tuple, T}) where {
        F, T <: Union{AbstractFloat, Complex}}
    return Base.FastMath.sqrt_fast(real(sum(abs2 ∘ f, u)))
end

@inline function NONLINEARSOLVE_DEFAULT_NORM(u::AbstractArray)
    return Base.FastMath.sqrt_fast(UNITLESS_ABS2(u))
end

@inline function NONLINEARSOLVE_DEFAULT_NORM(f::F, u::AbstractArray) where {F}
    return Base.FastMath.sqrt_fast(UNITLESS_ABS2(f, u))
end

@inline NONLINEARSOLVE_DEFAULT_NORM(u) = norm(u)
@inline NONLINEARSOLVE_DEFAULT_NORM(f::F, u) where {F} = norm(f.(u))
