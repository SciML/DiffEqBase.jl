using OrdinaryDiffEq, Test

f(u,p,t) = 2u
u0 = 0.5
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob,Tsit5())

function f(du, u, p, t)
    du .= 2.0 * u
end
prob = ODEProblem(f,u0,tspan)
@test_throws DiffEqBase.IncompatibleInitialConditionError sol = solve(prob,Tsit5())

prob = ODEProblem{false}(f,u0,tspan)
sol = solve(prob,Tsit5())

@test_throws DiffEqBase.ProblemSolverPairingError solve(prob, DFBDF())
@test_throws DiffEqBase.NoDefaultAlgorithmError solve(prob,nothing)
@test_throws DiffEqBase.NonSolverError solve(prob,5.0)

prob = ODEProblem{false}(f,u0,(nothing,nothing))
@test_throws DiffEqBase.NoTspanError solve(prob,Tsit5())