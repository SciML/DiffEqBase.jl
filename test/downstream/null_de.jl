using ModelingToolkit, OrdinaryDiffEq, SteadyStateDiffEq, Test

@variables t x(t) y(t)
eqs = [
    0 ~ x - y
    0 ~ y - x
]

@named sys = ODESystem(eqs, t)
sys = structural_simplify(sys)

prob = ODEProblem(sys, Pair[], (0, 10.0))
sol = solve(prob, Tsit5())

@test sol[x] == [0.0, 0.0]
@test sol[y] == [0.0, 0.0]

@variables t x y
eqs = [
    0 ~ x - y
    0 ~ y - x
]

@named sys = NonlinearSystem(eqs, [x, y], [])
sys = structural_simplify(sys)
prob = NonlinearProblem(sys,[])

sol = solve(prob, DynamicSS(Tsit5()))

@test sol[x] == 0.0
@test sol[y] == 0.0
