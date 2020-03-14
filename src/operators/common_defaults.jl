# The `update_coefficients!` interface
DEFAULT_UPDATE_FUNC(A,u,p,t) = A # no-op used by the basic operators
# isconstant(::AbstractDiffEqLinearOperator) = true # already defined in DiffEqBase
update_coefficients!(L::AbstractDiffEqLinearOperator,u,p,t) = L

# Routines that use the AbstractMatrix representation
Base.convert(::Type{AbstractArray}, L::AbstractDiffEqLinearOperator) = convert(AbstractMatrix, L)
Base.size(L::AbstractDiffEqLinearOperator, args...) = size(convert(AbstractMatrix,L), args...)
LinearAlgebra.opnorm(L::AbstractDiffEqLinearOperator, p::Real=2) = opnorm(convert(AbstractMatrix,L), p)
Base.getindex(L::AbstractDiffEqLinearOperator, i::Int) = convert(AbstractMatrix,L)[i]
Base.getindex(L::AbstractDiffEqLinearOperator, I::Vararg{Int, N}) where {N} =
  convert(AbstractMatrix,L)[I...]
for op in (:*, :/, :\)
  @eval Base.$op(L::AbstractDiffEqLinearOperator, x::Union{AbstractArray,Number}) = $op(convert(AbstractMatrix,L), x)
  @eval Base.$op(x::Union{AbstractVecOrMat,Number}, L::AbstractDiffEqLinearOperator) = $op(x, convert(AbstractMatrix,L))
end
LinearAlgebra.mul!(Y::AbstractArray, L::AbstractDiffEqLinearOperator, B::AbstractArray) =
  mul!(Y, convert(AbstractMatrix,L), B)
LinearAlgebra.ldiv!(Y::AbstractVecOrMat, L::AbstractDiffEqLinearOperator, B::AbstractVecOrMat) =
  ldiv!(Y, convert(AbstractMatrix,L), B)
for pred in (:isreal, :issymmetric, :ishermitian, :isposdef)
  @eval LinearAlgebra.$pred(L::AbstractDiffEqLinearOperator) = $pred(convert(AbstractArray, L))
end
for op in (:sum,:prod)
  @eval LinearAlgebra.$op(L::AbstractDiffEqLinearOperator; kwargs...) = $op(convert(AbstractArray, L); kwargs...)
end
LinearAlgebra.factorize(L::AbstractDiffEqLinearOperator) =
  FactorizedDiffEqArrayOperator(factorize(convert(AbstractMatrix, L)))
for fact in (:lu, :lu!, :qr, :qr!, :cholesky, :cholesky!, :ldlt, :ldlt!,
  :bunchkaufman, :bunchkaufman!, :lq, :lq!, :svd, :svd!)
  @eval LinearAlgebra.$fact(L::AbstractDiffEqLinearOperator, args...) =
    FactorizedDiffEqArrayOperator($fact(convert(AbstractMatrix, L), args...))
  @eval LinearAlgebra.$fact(L::AbstractDiffEqLinearOperator; kwargs...) =
    FactorizedDiffEqArrayOperator($fact(convert(AbstractMatrix, L); kwargs...))
end

# Routines that use the full matrix representation
Base.Matrix(L::AbstractDiffEqLinearOperator) = Matrix(convert(AbstractMatrix, L))
LinearAlgebra.exp(L::AbstractDiffEqLinearOperator) = exp(Matrix(L))
