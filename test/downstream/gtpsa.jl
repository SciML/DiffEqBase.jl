using OrdinaryDiffEq, ForwardDiff, GTPSA, Test

f!(du, u, p, t) = du .= p .* u

# Initial variables and parameters
x = [1.0, 2.0, 3.0]
p = [4.0, 5.0, 6.0]

prob = ODEProblem(f!, x, (0.0, 1.0), p)
sol = solve(prob, Tsit5(), reltol=1e-16, abstol=1e-16)

# Parametric GTPSA map
desc = Descriptor(3, 2, 3, 2) # 3 variables 3 parameters, both to 2nd order
dx = vars(desc)
dp = params(desc)
prob_GTPSA = ODEProblem(f!, x .+ dx, (0.0, 1.0), p .+ dp)
sol_GTPSA = solve(prob_GTPSA, Tsit5(), reltol=1e-16, abstol=1e-16)

@test sol.u[end] ≈ scalar.(sol_GTPSA.u[end]) # scalar gets 0th order part

# Compare Jacobian against ForwardDiff
J_FD = ForwardDiff.jacobian([x..., p...]) do t
    prob = ODEProblem(f!, t[1:3], (0.0, 1.0), t[4:6])
    sol = solve(prob, Tsit5(), reltol=1e-16, abstol=1e-16)
    sol.u[end]
end

@test J_FD ≈ GTPSA.jacobian(sol_GTPSA.u[end], include_params=true)

# Compare Hessians against ForwardDiff
for i in 1:3
    Hi_FD = ForwardDiff.hessian([x..., p...]) do t
        prob = ODEProblem(f!, t[1:3], (0.0, 1.0), t[4:6])
        sol = solve(prob, Tsit5(), reltol=1e-16, abstol=1e-16)
        sol.u[end][i]
    end
    @test Hi_FD ≈ GTPSA.hessian(sol_GTPSA.u[end][i], include_params=true)
end

