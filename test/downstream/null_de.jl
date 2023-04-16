using ModelingToolkit, OrdinaryDiffEq, SteadyStateDiffEq, Test

@variables t x(t) y(t)
eqs = [0 ~ x - y
       0 ~ y - x]

@named sys = ODESystem(eqs, t)
sys = structural_simplify(sys)

prob = ODEProblem(sys, Pair[], (0, 10.0))
sol = solve(prob, Tsit5())

@test sol[x] == [0.0, 0.0]
@test sol[y] == [0.0, 0.0]

for kwargs in [
    Dict(:saveat => 0:0.1:1),
    Dict(:save_start => false),
    Dict(:save_end => false),
]
    sol = solve(prob, kwargs...)
    init_integ = init(prob, kwargs...)
    solve!(init_integ)
    step_integ = init(prob, kwargs...)
    step!(step_integ, prob.tspan[end] - prob.tspan[begin])
    @test sol.u[end] == init_integ.u
    @test sol.t[end] == init_integ.t
    @test sol.u[end] == step_integ.u
    @test sol.t[end] == step_integ.t
end

@variables t x y
eqs = [0 ~ x - y
       0 ~ y - x]

@named sys = NonlinearSystem(eqs, [x, y], [])
sys = structural_simplify(sys)
prob = NonlinearProblem(sys, [])

sol = solve(prob, DynamicSS(Tsit5()))

@test sol[x] == 0.0
@test sol[y] == 0.0
