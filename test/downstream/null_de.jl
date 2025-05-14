using ModelingToolkit, OrdinaryDiffEq, SteadyStateDiffEq, Test
using StaticArrays
using ModelingToolkit: t_nounits as t, D_nounits as D
using ForwardDiff

@variables x(t) y(t)
eqs = [0 ~ x - y
       0 ~ y - x]

@named sys = ODESystem(eqs, t)
sys = structural_simplify(sys)

prob = ODEProblem(sys, Pair[], (0, 10.0))
sol = solve(prob, Tsit5())

@test sol.dense
@test sol[x] == [0.0, 0.0]
@test sol[y] == [0.0, 0.0]

for kwargs in [
    Dict(:saveat => 0:0.1:1),
    Dict(:save_start => false),
    Dict(:save_end => false)
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

@variables x y
eqs = [0 ~ x - y
       0 ~ y - x]

@named sys = NonlinearSystem(eqs, [x, y], [])
sys = structural_simplify(sys)
prob = NonlinearProblem(sys, [])

sol = solve(prob, DynamicSS(Tsit5()))

@test sol[x] == 0.0
@test sol[y] == 0.0

# https://github.com/SciML/NonlinearSolve.jl/issues/387

using NonlinearSolve
function unsat(du, u, p)
    du[1] = 1
end
unsat_f = NonlinearFunction(unsat; resid_prototype = zeros(1))
unsatprob = NonlinearLeastSquaresProblem(unsat_f, nothing)
sol = solve(unsatprob) # Success
@test sol.retcode == SciMLBase.ReturnCode.Failure
@test sol.resid == [1.0]

# Issue#2664
@testset "remake type promotion with empty initial conditions" begin
    @parameters P
    @variables x(t)
    # numerical ODE: x′(t) = P with x(0) = 0
    sys_num = structural_simplify(ODESystem([D(x) ~ P], t, [x], [P]; name = :sys))
    prob_num_uninit = ODEProblem(sys_num, [x => 0.0], (0.0, 1.0), [P => NaN]) # uninitialized problem
    x_at_1_num(P) = solve(remake(prob_num_uninit; p = [sys_num.P => P]), Tsit5())(
        1.0, idxs = x)

    # analytical solution: x(t) = P*t
    sys_anal = structural_simplify(ODESystem([x ~ P * t], t, [x], [P]; name = :sys))
    prob_anal_uninit = ODEProblem(sys_anal, [], (0.0, 1.0), [P => NaN])
    x_at_1_anal(P) = solve(remake(prob_anal_uninit; p = [sys_anal.P => P]), Tsit5())(
        1.0, idxs = x)

    @test_nowarn x_at_1_num(1.0)
    @test_nowarn x_at_1_anal(1.0)
    @test_nowarn ForwardDiff.derivative(x_at_1_num, 1.0)
    @test_nowarn ForwardDiff.derivative(x_at_1_anal, 1.0)
end

@testset "Null OOP NonlinearLeastSquaresProblem" begin
    fn = NonlinearFunction{false}(; resid_prototype = SA[1.0]) do u, p
        return zero(u)
    end
    prob = NonlinearLeastSquaresProblem(fn, nothing)
    sol = solve(prob)
    @test sol.resid isa SVector{1, Float64}
end
