### AbstractDiffEqOperator interface

#=
1. Function call and multiplication: A(t,u,du) for inplace and du = A(t,u) for
   out-of-place, meaning A*u and A_mul_B!.
2. If the operator is not a constant, update it with (t,u). A mutating form, i.e.
   update_coefficients!(A,t,u) that changes the internal coefficients, and a
   out-of-place form B = update_coefficients(A,t,u).
3. is_constant(A) trait for whether the operator is constant or not.
4. diagonal, symmetric, etc traits from LinearMaps.jl.
5. Optional: expm(A). Required for simple exponential integration.
6. Optional: expmv(A,t,u) = expm(t*A)*u and expmv!(v,A::DiffEqOperator,t,u)
   Required for sparse-saving exponential integration.
7. Optional: factorizations. A_ldiv_B, factorize et. al. This is only required
   for algorithms which use the factorization of the operator (Crank-Nicholson),
   and only for when the default linear solve is used.
=#

# Inherits the standard assumptions of an AbstractLinearMap
# Extra standard assumptions
is_constant(A::AbstractDiffEqOperator) = true
update_coefficients!(A,t,u) = nothing
update_coefficients(A,t,u) = A

# Generic fallbacks
Base.expm(A::AbstractDiffEqOperator,t) = expm(t*A)
expmv(A::AbstractDiffEqOperator,t,u) = expm(t,A)*u
expmv!(v,A::AbstractDiffEqOperator,t,u) = A_mul_B!(v,expm(t,A),u)

### Constant DiffEqOperator defined by an array
struct DiffEqArrayOperator{T,Arr<:AbstractMatrix{T},F} <: AbstractDiffEqOperator{T}
  A::Arr
  _isreal::Bool
  _issymmetric::Bool
  _ishermitian::Bool
  _isposdef::Bool
  update_func::F
end
DEFAULT_UPDATE_FUNC = (A,t,u)->nothing
DiffEqArrayOperator{T}(A::AbstractMatrix{T},update_func = DEFAULT_UPDATE_FUNC) = DiffEqArrayOperator{T,typeof(A)}(
                       A,isreal(A),issymmetric(A),ishermitian(A),isposdef(A),update_func)

Base.isreal(A::DiffEqArrayOperator) = A._isreal
Base.issymmetric(A::DiffEqArrayOperator) = A._issymmetric
Base.ishermitian(A::DiffEqArrayOperator) = A._ishermitian
Base.isposdef(A::DiffEqArrayOperator) = A._isposdef

update_coefficients!(A::DiffEqArrayOperator,t,u) = A.update_func(A.A,t,u)
update_coefficients(A::DiffEqArrayOperator,t,u) = (A.update_func(A.A,t,u); A.A)

function (A::DiffEqArrayOperator)(t,u)
  update_coefficients!(A,t,u)
  A*u
end
function (A::DiffEqArrayOperator)(t,u,du)
  update_coefficients!(A,t,u)
  A_mul_B!(du,A.A,u)
end

### Forward some extra operations
Base.:*(A::DiffEqArrayOperator,B::DiffEqArrayOperator) = A.A*B.A
Base.:*(A::DiffEqArrayOperator,B::Number) = A.A*B
Base.:*(A::DiffEqArrayOperator,B::AbstractVector) = A.A*B
Base.:*(A::DiffEqArrayOperator,B::AbstractMatrix) = A.A*B
Base.:*(A,B::DiffEqArrayOperator) = A*B.A
Base.A_mul_B!(v,A::DiffEqArrayOperator,B::DiffEqArrayOperator) = A_mul_B!(v,A.A,B.A)
Base.A_mul_B!(v,A::DiffEqArrayOperator,B::Number) = A_mul_B!(v,A.A,B)
Base.A_mul_B!(v,A::DiffEqArrayOperator,B::AbstractVector) = A_mul_B!(v,A.A,B)
Base.A_mul_B!(v,A::DiffEqArrayOperator,B::AbstractMatrix) = A_mul_B!(v,A.A,B)
Base.A_mul_B!(v,A,B::DiffEqArrayOperator) = A_mul_B!(v,A,B.A)
Base.:+(A::DiffEqArrayOperator,B::DiffEqArrayOperator) = A.A+B.A
Base.:+(A::DiffEqArrayOperator,B) = A.A+B
Base.:+(B,A::DiffEqArrayOperator) = B+A.A
Base.:-(A::DiffEqArrayOperator,B::DiffEqArrayOperator) = A.A-B.A
Base.:-(A::DiffEqArrayOperator,B) = A.A-B
Base.:-(B,A::DiffEqArrayOperator) = B-A.A
Base.:/(A::DiffEqArrayOperator,B::DiffEqArrayOperator) = A.A/B.A
Base.:/(A::DiffEqArrayOperator,B) = A.A/B
Base.:/(B,A::DiffEqArrayOperator) = B/A.A
Base.:\(A::DiffEqArrayOperator,B::DiffEqArrayOperator) = A.A\B.A
Base.:\(A::DiffEqArrayOperator,B) = A.A\B
Base.:\(B,A::DiffEqArrayOperator) = B\A.A
Base.expm(A::DiffEqArrayOperator) = expm(A.A)
Base.A_ldiv_B!(Y,A::DiffEqArrayOperator, B::DiffEqArrayOperator) = A_ldiv_B!(Y,A.A, B.B)
Base.A_ldiv_B!(Y,A::DiffEqArrayOperator, B) = A_ldiv_B!(Y,A.A, B)
Base.A_ldiv_B!(Y,A, B::DiffEqArrayOperator) = A_ldiv_B!(Y,A, B.A)
Base.factorize(A::DiffEqArrayOperator) = factorize(A.A)
Base.lufact(A::DiffEqArrayOperator,args...)    = lufact(A.A,args...)
Base.lufact!(A::DiffEqArrayOperator,args...)   = lufact!(A.A,args...)
Base.qrfact(A::DiffEqArrayOperator,args...)    = qrfact(A.A,args...)
Base.qrfact!(A::DiffEqArrayOperator,args...)   = qrfact!(A.A,args...)
Base.cholfact(A::DiffEqArrayOperator,args...)  = cholfact(A.A,args...)
Base.cholfact!(A::DiffEqArrayOperator,args...) = cholfact!(A.A,args...)
Base.ldltfact(A::DiffEqArrayOperator,args...)  = ldltfact(A.A,args...)
Base.ldltfact!(A::DiffEqArrayOperator,args...) = ldltfact!(A.A,args...)
Base.bkfact(A::DiffEqArrayOperator,args...)    = bkfact(A.A,args...)
Base.bkfact!(A::DiffEqArrayOperator,args...)   = bkfact!(A.A,args...)
Base.lqfact(A::DiffEqArrayOperator,args...)    = lqfact(A.A,args...)
Base.lqfact!(A::DiffEqArrayOperator,args...)   = lqfact!(A.A,args...)
Base.svdfact(A::DiffEqArrayOperator,args...)   = svdfact(A.A,args...)
Base.svdfact!(A::DiffEqArrayOperator,args...)  = svdfact!(A.A,args...)
Base.size(A::DiffEqArrayOperator) = size(A.A)
@inline Base.getindex(A::DiffEqArrayOperator,i::Int) = A.A[i]
@inline Base.getindex{N}(A::DiffEqArrayOperator,I::Vararg{Int, N}) = A.A[I]
@inline Base.setindex!(A::DiffEqArrayOperator, v, i::Int) = (A.A[i]=v)
@inline Base.setindex!{N}(A::DiffEqArrayOperator, v, I::Vararg{Int, N}) = (A.A[I]=v)
