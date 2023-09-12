using DiffEqBase
using DiffEqBase: InternalFalsi, InternalITP, IntervalNonlinearProblem
using ForwardDiff

for Rootfinder in (InternalFalsi, InternalITP)
    rf = Rootfinder()
    # From SimpleNonlinearSolve
    f = (u, p) -> u * u - p
    tspan = (1.0, 20.0)
    g = function (p)
        probN = IntervalNonlinearProblem{false}(f, typeof(p).(tspan), p)
        sol = solve(probN, rf)
        return sol.u
    end

    for p in (1.0,) #1.1:0.1:100.0
        @test g(p) ≈ sqrt(p)
        #@test ForwardDiff.derivative(g, p) ≈ 1 / (2 * sqrt(p))
    end

    # https://github.com/SciML/DiffEqBase.jl/issues/916 
    inp = IntervalNonlinearProblem((t, p) -> min(-1.0 + 0.001427344607477125 * t, 1e-9),
        (699.0079267259368, 700.6176418816023))
    @test solve(inp, rf).u ≈ 700.6016590257979

    # Flipped signs & reversed tspan test for bracketing algorithms
    f1(u, p) = u * u - p
    f2(u, p) = p - u * u
    for p in 1:4
        inp1 = IntervalNonlinearProblem(f1, (1.0, 2.0), p)
        inp2 = IntervalNonlinearProblem(f2, (1.0, 2.0), p)
        inp3 = IntervalNonlinearProblem(f1, (2.0, 1.0), p)
        inp4 = IntervalNonlinearProblem(f2, (2.0, 1.0), p)
        @test abs.(solve(inp1, rf).u) ≈ sqrt.(p)
        @test abs.(solve(inp2, rf).u) ≈ sqrt.(p)
        @test abs.(solve(inp3, rf).u) ≈ sqrt.(p)
        @test abs.(solve(inp4, rf).u) ≈ sqrt.(p)
    end
end
