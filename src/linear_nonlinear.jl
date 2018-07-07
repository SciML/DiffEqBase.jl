using LinearAlgebra

mutable struct LinSolveFactorize{F}
  factorization::F
  A
end
LinSolveFactorize(factorization) = LinSolveFactorize(factorization,nothing)
function (p::LinSolveFactorize)(x,A,b,update_matrix=false)
  if update_matrix
    p.A = p.factorization(A)
  end
  ldiv!(x,p.A,b)
end
DEFAULT_LINSOLVE = LinSolveFactorize(lu!)
function (p::LinSolveFactorize)(::Type{Val{:init}},f,u0_prototype)
  LinSolveFactorize(p.factorization,nothing)
end
