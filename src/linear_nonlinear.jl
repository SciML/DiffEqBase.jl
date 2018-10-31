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
DEFAULT_LINSOLVE = LinSolveFactorize(lu!)
function (p::LinSolveFactorize)(::Type{Val{:init}},f,u0_prototype)
  LinSolveFactorize(p.factorization,nothing)
end
