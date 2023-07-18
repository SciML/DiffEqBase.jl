using DiffEqBase
using DiffEqBase: InternalFalsi, InternalITP, IntervalNonlinearProblem
using ForwardDiff

for Rootfinder in (InternalFalsi, InternalITP)
    # From SimpleNonlinearSolve
    f = (u, p) -> u * u - p
    tspan = (1.0, 20.0)
    g = function (p)
        probN = IntervalNonlinearProblem{false}(f, typeof(p).(tspan), p)
        sol = solve(probN, Rootfinder())
        return sol.u
    end

    for p in (1.0,) #1.1:0.1:100.0
        @test g(p) ≈ sqrt(p)
        #@test ForwardDiff.derivative(g, p) ≈ 1 / (2 * sqrt(p))
    end

    # https://github.com/SciML/DiffEqBase.jl/issues/916 
    inp = IntervalNonlinearProblem((t, p) -> min(-1.0 + 0.001427344607477125 * t, 1e-9),
                                   (699.0079267259368, 700.6176418816023))
    @test solve(inp, Rootfinder()).u ≈ 700.6016590257979
end