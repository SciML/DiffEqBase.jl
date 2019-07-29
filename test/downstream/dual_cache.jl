using LinearAlgebra, OrdinaryDiffEq
function foo(du, u, (A, tmp), t)
    tmp = DiffEqBase.get_tmp(tmp, u)
    mul!(tmp, A, u)
    @. du = u + tmp
    nothing
end
prob = ODEProblem(foo, ones(5, 5), (0., 1.0), (ones(5,5), DiffEqBase.dualcache(zeros(5,5))))
solve(prob, TRBDF2())
