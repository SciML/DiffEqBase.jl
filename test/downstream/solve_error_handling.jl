using OrdinaryDiffEq, StochasticDiffEq, Test, Sundials

f(u, p, t) = 2u
u0 = 0.5
tspan = (0.0, 1.0)
prob = ODEProblem(f, u0, tspan)
sol = solve(prob, Tsit5())

function f(du, u, p, t)
    du .= 2.0 * u
end
prob = ODEProblem(f, u0, tspan)
@test_throws DiffEqBase.IncompatibleInitialConditionError sol=solve(prob, Tsit5())

prob = ODEProblem{false}(f, u0, tspan)
sol = solve(prob, Tsit5())
sol = solve(prob, nothing, alg = Tsit5())
sol = init(prob, nothing, alg = Tsit5())

prob = ODEProblem{false}(f, 1.0 + im, tspan)
@test_throws DiffEqBase.ComplexSupportError solve(prob, CVODE_Adams())

@test_throws DiffEqBase.ProblemSolverPairingError solve(prob, DFBDF())
@test_throws DiffEqBase.NoDefaultAlgorithmError solve(prob, nothing)
@test_throws DiffEqBase.NoDefaultAlgorithmError init(prob, nothing)
@test_throws DiffEqBase.NonSolverError solve(prob, 5.0)

prob = ODEProblem{false}(f, u0, (nothing, nothing))
@test_throws DiffEqBase.NoTspanError solve(prob, Tsit5())

prob = ODEProblem{false}(f, u0, (nothing, nothing))
@test_throws DiffEqBase.NaNTspanError solve(prob, Tsit5())

prob = ODEProblem{false}(f, Any[1.0, 1.0f0], tspan)
@test_throws DiffEqBase.NonConcreteEltypeError solve(prob, Tsit5())

prob = ODEProblem{false}(f, (1.0, 1.0f0), tspan)
@test_throws DiffEqBase.TupleStateError solve(prob, Tsit5())

prob = ODEProblem{false}(f, u0, (0.0 + im, 1.0))
@test_throws DiffEqBase.ComplexTspanError solve(prob, Tsit5())

for u0 in ([0.0, 0.0], nothing)
    fmm = ODEFunction(f, mass_matrix = zeros(3, 3))
    prob = ODEProblem(fmm, u0, (0.0, 1.0))
    @test_throws DiffEqBase.IncompatibleMassMatrixError solve(prob, Tsit5())
end

# Allow empty mass matrix for empty u0
fmm = ODEFunction((du, u, t) -> nothing, mass_matrix = zeros(0, 0))
prob = ODEProblem(fmm, nothing, (0.0, 1.0))
sol = solve(prob, Tsit5())
@test isa(sol, DiffEqBase.ODESolution)

f(du, u, p, t) = du .= 1.01u
function g(du, u, p, t)
    du[1, 1] = 0.3u[1]
    du[1, 2] = 0.6u[1]
    du[1, 3] = 0.9u[1]
    du[1, 4] = 0.12u[1]
    du[2, 1] = 1.2u[2]
    du[2, 2] = 0.2u[2]
    du[2, 3] = 0.3u[2]
    du[2, 4] = 1.8u[2]
end

prob = SDEProblem(f,
    g,
    randn(ComplexF64, 2),
    (0.0, 1.0),
    noise_rate_prototype = complex(zeros(2, 4)),
    noise = StochasticDiffEq.RealWienerProcess(0.0, zeros(3)))
@test_throws DiffEqBase.NoiseSizeIncompatabilityError solve(prob, LambaEM())
