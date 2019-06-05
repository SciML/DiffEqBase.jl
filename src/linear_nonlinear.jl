using LinearAlgebra, SparseArrays, SuiteSparse
mutable struct LinSolveFactorize{F}
  factorization::F
  A
end
LinSolveFactorize(factorization) = LinSolveFactorize(factorization,nothing)
function (p::LinSolveFactorize)(x,A,b,update_matrix=false)
  if update_matrix
    p.A = p.factorization(A)
  end
  if typeof(p.A) <: SuiteSparse.UMFPACK.UmfpackLU || typeof(p.factorization) <: typeof(lu)
    ldiv!(x,p.A,b) # No 2-arg form for SparseArrays!
  else
    x .= b
    ldiv!(p.A,x)
  end
end
function (p::LinSolveFactorize)(::Type{Val{:init}},f,u0_prototype)
  LinSolveFactorize(p.factorization,nothing)
end

### Default Linsolve

# Try to be as smart as possible
# lu! if Matrix
# lu if sparse
# gmres if operator

mutable struct DefaultLinSolve
  A
end
DefaultLinSolve() = DefaultLinSolve(nothing)

function (p::DefaultLinSolve)(x,A,b,update_matrix=false;kwargs...)
  if update_matrix
    if typeof(A) <: Matrix
      blasvendor = BLAS.vendor()
      if (blasvendor === :openblas || blasvendor === :openblas64) && size(A,1) <= 500 # if the user doesn't use OpenBLAS, we assume that is a much better BLAS implementation like MKL
        p.A = RecursiveFactorization.lu!(A)
      else
        p.A = lu!(A)
      end
    elseif typeof(A) <: SparseMatrixCSC
      p.A = factorize(A)
    elseif !(typeof(A) <: AbstractDiffEqOperator)
      # Most likely QR is the one that is overloaded
      # Works on things like CuArrays
      p.A = qr(A)
    end
  end

  if typeof(A) <: Matrix # No 2-arg form for SparseArrays!
    x .= b
    ldiv!(p.A,x)
  # Missing a little bit of efficiency in a rare case
  #elseif typeof(A) <: DiffEqArrayOperator
  #  ldiv!(x,p.A,b)
  elseif typeof(A) <: AbstractDiffEqOperator
    # No good starting guess, so guess zero
    x .= false
    IterativeSolvers.gmres!(x,A,b,initially_zero=true)
  else
    ldiv!(x,p.A,b)
  end
end

function (p::DefaultLinSolve)(::Type{Val{:init}},f,u0_prototype)
  DefaultLinSolve()
end

const DEFAULT_LINSOLVE = DefaultLinSolve()

### Default GMRES

# Easily change to GMRES

mutable struct LinSolveGMRES{A}
  iterable
  kwargs::A
end
LinSolveGMRES(;kwargs...) = LinSolveGMRES(nothing, kwargs)

function (f::LinSolveGMRES)(x,A,b,update_matrix=false; tol, kwargs...)
  if f.iterable === nothing
    f.iterable = IterativeSolvers.gmres_iterable!(x,A,b;initially_zero=true,restart=5,maxiter=5,tol=1e-16,f.kwargs...,kwargs...)
  end
  x .= false
  iter = f.iterable
  for residual in iter
    residual ≤ tol && break # only use absolute tolerance
  end
  return nothing
end

function (p::LinSolveGMRES)(::Type{Val{:init}},f,u0_prototype)
  LinSolveGMRES(nothing, p.kwargs)
end

# scaling for iterative solvers
struct ScaleVector{T}
  x::T
  isleft::Bool
end
function LinearAlgebra.ldiv!(v::ScaleVector, x)
  if v.isleft
    return @.. x = x * v.x
  else
    return @.. x = x / v.x
  end
end
function LinearAlgebra.ldiv!(y, v::ScaleVector, x)
  if v.isleft
    return @.. y = x * v.x
  else
    return @.. y = x / v.x
  end
end
