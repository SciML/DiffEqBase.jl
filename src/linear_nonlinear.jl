type LinSolveFactorize{F}
  factorization::F
  A
end
LinSolveFactorize(factorization) = LinSolveFactorize(factorization,nothing)
function (p::LinSolveFactorize)(x,A,b,update_matrix=false)
  if update_matrix
    p.A = p.factorization(A)
  end
  A_ldiv_B!(x,p.A,b)
end
DEFAULT_LINSOLVE = LinSolveFactorize(lufact!)
