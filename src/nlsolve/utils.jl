## Anderson acceleration

"""
    anderson(gz, cache, integrator)

Return the value for the next fixed-point iteration by performing Anderson acceleration
based on the current evaluation `gz = G(z)` of the fixed-point iteration, and the settings
and history in the `cache`.
"""
@muladd function anderson(gz, cache, integrator)
  @unpack dz,Δgzs,Q,R,γs,history,droptol = cache

  # increase size of history
  history += 1
  
  # remove oldest history if maximum size is exceeded
  max_history = length(Δgzs)
  if history > max_history
    # circularly shift differences of G(z)
    for i in 1:(max_history-1)
      Δgzs[i] = Δgzs[i + 1]
    end
  
    # delete left-most column of QR decomposition
    qrdelete!(Q, R, max_history)
  
    # update size of history
    history = max_history
  end
  
  # update history of differences of G(z)
  Δgzs[history] = @.. gz - cache.gzprev
  
  # replace/add difference of residuals as right-most column to QR decomposition
  qradd!(Q, R, _vec(dz .- cache.dzprev), history)
  
  # update cached values
  cache.dzprev = dz
  cache.gzprev = gz
  
  # define current Q and R matrices
  Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))
  
  # check condition (TODO: incremental estimation)
  if droptol !== nothing
    while cond(R) > droptol && history > 1
      qrdelete!(Q, R, history)
      history -= 1
      Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))
    end
  end
  
  # solve least squares problem
  γscur = view(γs, 1:history)
  ldiv!(Rcur, mul!(γscur, Qcur', _vec(dz)))
  if has_destats(integrator)
    integrator.destats.nsolve += 1
  end
  
  # update iterate
  for i in 1:history
    gz = @.. gz - γs[i] * Δgzs[i]
  end
  
  # update cached values
  cache.history = history

  gz
end

"""
    anderson!(gz, cache, integrator)

Update the current evaluation `gz = G(z)` of the fixed-point iteration in-place by
performing Anderson acceleration based on the settings and history in the `cache`.
"""
@muladd function anderson!(gz, cache, integrator)
  @unpack gzprev,dz,dzprev,Δgzs,Q,R,γs,history,droptol = cache

  # increase size of history
  history += 1

  # remove oldest history if maximum size is exceeded
  max_history = length(Δgzs)
  if history > max_history
    # circularly shift differences of z
    ptr = Δgzs[1]
    for i in 1:(max_history-1)
      Δgzs[i] = Δgzs[i + 1]
    end
    Δgzs[max_history] = ptr

    # delete left-most column of QR decomposition
    qrdelete!(Q, R, max_history)

    # update size of history
    history = max_history
  end

  # update history of differences of g(z)
  @.. Δgzs[history] = gz - gzprev

  # replace/add difference of residuals as right-most column to QR decomposition
  @.. dzprev = dz - dzprev
  qradd!(Q, R, vec(dzprev), history)

  # update cached values
  @.. dzprev = dz
  @.. gzprev = gz

  # define current Q and R matrices
  Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))

  # check condition (TODO: incremental estimation)
  if droptol !== nothing
    while cond(R) > droptol && history > 1
      qrdelete!(Q, R, history)
      history -= 1
      Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))
    end
  end

  # solve least squares problem
  γscur = view(γs, 1:history)
  ldiv!(Rcur, mul!(γscur, Qcur', vec(dz)))
  if has_destats(integrator)
    integrator.destats.nsolve += 1
  end

  # update next iterate
  for i in 1:history
    @.. gz = gz - γs[i] * Δgzs[i]
  end

  # update cached values
  cache.history = history

  nothing
end

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