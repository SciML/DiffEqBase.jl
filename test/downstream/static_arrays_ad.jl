using OrdinaryDiffEq, ForwardDiff, StaticArrays, Test
f(u, p, t) = copy(u)

du1 = ForwardDiff.derivative(5.0) do x
    prob = ODEProblem(f, [x], (0.0, 1.0), nothing)
    sol = solve(prob, Tsit5())
    sol.u[end]
end

du2 = ForwardDiff.derivative(5.0) do x
    prob = ODEProblem(f, SVector(x), (0.0, 1.0), nothing)
    sol = solve(prob, Tsit5())
    sol.u[end]
end

@test du1 ≈ du2
