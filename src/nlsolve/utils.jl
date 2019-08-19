## helpers for Anderson acceleration

"""
    qrdelete!(Q, R, k)

Delete the left-most column of F = Q[:, 1:k] * R[1:k, 1:k] by updating Q and R.
Only Q[:, 1:(k-1)] and R[1:(k-1), 1:(k-1)] are valid on exit.
"""
function qrdelete!(Q::AbstractMatrix, R::AbstractMatrix, k::Int)
  n, m = size(Q)
  m == LinearAlgebra.checksquare(R) || throw(DimensionMismatch())
  1 ≤ k ≤ m || throw(ArgumentError())

  # apply Givens rotations
  for i in 2:k
      g = first(givens(R, i - 1, i, i))
      lmul!(g, R)
      rmul!(Q, g')
  end

  # move columns of R
  @inbounds for j in 1:(k-1)
    for i in 1:(k-1)
      R[i, j] = R[i, j + 1]
    end
  end

  Q, R
end

"""
    qradd!(Q, R, v, k)

Replace the right-most column of F = Q[:, 1:k] * R[1:k, 1:k] with v by updating Q and R.
This implementation modifies vector v as well. Only Q[:, 1:k] and R[1:k, 1:k] are valid on
exit.
"""
function qradd!(Q::AbstractMatrix, R::AbstractMatrix, v::AbstractVector, k::Int)
  n, m = size(Q)
  n == length(v) || throw(DimensionMismatch())
  m == LinearAlgebra.checksquare(R) || throw(DimensionMismatch())
  1 ≤ k ≤ m || throw(ArgumentError())

  @inbounds for i in 1:(k-1)
    q = view(Q, :, i)
    r = dot(q, v)

    R[i, k] = r
    axpy!(-r, q, v)
  end

  @inbounds begin
    d = norm(v)
    R[k, k] = d
    @.. @view(Q[:, k]) = v / d
  end

  Q, R
end

function qradd!(Q::AbstractMatrix, R::AbstractMatrix, v::Number, k::Int)
  1 == LinearAlgebra.checksquare(Q) == LinearAlgebra.checksquare(R) ||
    throw(DimensionMismatch())
  k == 1 || throw(ArgumentError())

  R[1, 1] = abs(v)
  Q[1, 1] = one(v)

  Q, R
end