"""
$(TYPEDEF)

TODO
"""
struct DiffEqIdentity{T,N} <: AbstractDiffEqLinearOperator{T} end

DiffEqIdentity(u) = DiffEqIdentity{eltype(u),length(u)}()
Base.size(::DiffEqIdentity{T,N}) where {T,N} = (N,N)
Base.size(::DiffEqIdentity{T,N}, m::Integer) where {T,N} = (m == 1 || m == 2) ? N : 1
LinearAlgebra.opnorm(::DiffEqIdentity{T,N}, p::Real=2) where {T,N} = one(T)
Base.convert(::Type{AbstractMatrix}, ::DiffEqIdentity{T,N}) where {T,N} =
                                              LinearAlgebra.Diagonal(ones(T,N))

for op in (:*, :/, :\)
  @eval Base.$op(::DiffEqIdentity{T,N}, x::AbstractVecOrMat) where {T,N} = $op(I, x)
  @eval Base.$op(::DiffEqIdentity{T,N}, x::AbstractArray) where {T,N} = $op(I, x)
  @eval Base.$op(x::AbstractVecOrMat, ::DiffEqIdentity{T,N}) where {T,N} = $op(x, I)
  @eval Base.$op(x::AbstractArray, ::DiffEqIdentity{T,N}) where {T,N} = $op(x, I)

end
LinearAlgebra.mul!(Y::AbstractVecOrMat, ::DiffEqIdentity, B::AbstractVecOrMat) = Y .= B
LinearAlgebra.ldiv!(Y::AbstractVecOrMat, ::DiffEqIdentity, B::AbstractVecOrMat) = Y .= B

LinearAlgebra.mul!(Y::AbstractArray, ::DiffEqIdentity, B::AbstractArray) = Y .= B
LinearAlgebra.ldiv!(Y::AbstractArray, ::DiffEqIdentity, B::AbstractArray) = Y .= B

for pred in (:isreal, :issymmetric, :ishermitian, :isposdef)
  @eval LinearAlgebra.$pred(::DiffEqIdentity) = true
end

"""
    DiffEqScalar(val[; update_func])

Represents a time-dependent scalar/scaling operator. The update function
is called by `update_coefficients!` and is assumed to have the following
signature:

    update_func(oldval,u,p,t) -> newval

You can also use `setval!(α,val)` to bypass the `update_coefficients!`
interface and directly mutate the scalar's value.
"""
mutable struct DiffEqScalar{T<:Number,F} <: AbstractDiffEqLinearOperator{T}
  val::T
  update_func::F
  DiffEqScalar(val::T; update_func=DEFAULT_UPDATE_FUNC) where {T} =
    new{T,typeof(update_func)}(val, update_func)
end

Base.convert(::Type{Number}, α::DiffEqScalar) = α.val
Base.convert(::Type{DiffEqScalar}, α::Number) = DiffEqScalar(α)
Base.size(::DiffEqScalar) = ()
Base.size(::DiffEqScalar, ::Integer) = 1
update_coefficients!(α::DiffEqScalar,u,p,t) = (α.val = α.update_func(α.val,u,p,t); α)
setval!(α::DiffEqScalar, val) = (α.val = val; α)
isconstant(α::DiffEqScalar) = α.update_func == DEFAULT_UPDATE_FUNC

for op in (:*, :/, :\)
  @eval Base.$op(α::DiffEqScalar, x::Union{AbstractArray,Number}) = $op(α.val, x)
  @eval Base.$op(x::Union{AbstractArray,Number}, α::DiffEqScalar) = $op(x, α.val)
  @eval Base.$op(x::DiffEqScalar, y::DiffEqScalar) = $op(x.val, y.val)
end

for op in (:-, :+)
  @eval Base.$op(α::DiffEqScalar, x::Number) = $op(α.val, x)
  @eval Base.$op(x::Number, α::DiffEqScalar) = $op(x, α.val)
  @eval Base.$op(x::DiffEqScalar, y::DiffEqScalar) = $op(x.val, y.val)
end

LinearAlgebra.lmul!(α::DiffEqScalar, B::AbstractArray) = lmul!(α.val, B)
LinearAlgebra.rmul!(B::AbstractArray, α::DiffEqScalar) = rmul!(B, α.val)
LinearAlgebra.mul!(Y::AbstractArray, α::DiffEqScalar, B::AbstractArray) = mul!(Y, α.val, B)
LinearAlgebra.axpy!(α::DiffEqScalar, X::AbstractArray, Y::AbstractArray) = axpy!(α.val, X, Y)
Base.abs(α::DiffEqScalar) = abs(α.val)

"""
    DiffEqArrayOperator(A[; update_func])

Represents a time-dependent linear operator given by an AbstractMatrix. The
update function is called by `update_coefficients!` and is assumed to have
the following signature:

    update_func(A::AbstractMatrix,u,p,t) -> [modifies A]

You can also use `setval!(α,A)` to bypass the `update_coefficients!` interface
and directly mutate the array's value.
"""
mutable struct DiffEqArrayOperator{T,AType<:AbstractMatrix{T},F} <: AbstractDiffEqLinearOperator{T}
  A::AType
  update_func::F
  DiffEqArrayOperator(A::AType; update_func=DEFAULT_UPDATE_FUNC) where {AType} =
    new{eltype(A),AType,typeof(update_func)}(A, update_func)
end

update_coefficients!(L::DiffEqArrayOperator,u,p,t) = (L.update_func(L.A,u,p,t); L)
setval!(L::DiffEqArrayOperator, A) = (L.A = A; L)
isconstant(L::DiffEqArrayOperator) = L.update_func == DEFAULT_UPDATE_FUNC

Base.convert(::Type{AbstractMatrix}, L::DiffEqArrayOperator) = L.A
Base.setindex!(L::DiffEqArrayOperator, v, i::Int) = (L.A[i] = v)
Base.setindex!(L::DiffEqArrayOperator, v, I::Vararg{Int, N}) where {N} = (L.A[I...] = v)

"""
    FactorizedDiffEqArrayOperator(F)

Like DiffEqArrayOperator, but stores a Factorization instead.

Supports left division and `ldiv!` when applied to an array.
"""
struct FactorizedDiffEqArrayOperator{T<:Number,FType<:Factorization{T}} <: AbstractDiffEqLinearOperator{T}
  F::FType
end

Base.Matrix(L::FactorizedDiffEqArrayOperator) = Matrix(L.F)
Base.convert(::Type{AbstractMatrix}, L::FactorizedDiffEqArrayOperator) = convert(AbstractMatrix, L.F)
Base.size(L::FactorizedDiffEqArrayOperator, args...) = size(L.F, args...)
LinearAlgebra.ldiv!(Y::AbstractVecOrMat, L::FactorizedDiffEqArrayOperator, B::AbstractVecOrMat) = ldiv!(Y, L.F, B)
LinearAlgebra.ldiv!(L::FactorizedDiffEqArrayOperator, B::AbstractVecOrMat) = ldiv!(L.F, B)
Base.:\(L::FactorizedDiffEqArrayOperator, x::AbstractVecOrMat) = L.F \ x

# The (u,p,t) and (du,u,p,t) interface
for T in [DiffEqScalar, DiffEqArrayOperator, FactorizedDiffEqArrayOperator, DiffEqIdentity]
  (L::T)(u,p,t) = (update_coefficients!(L,u,p,t); L * u)
  (L::T)(du,u,p,t) = (update_coefficients!(L,u,p,t); mul!(du,L,u))
end
