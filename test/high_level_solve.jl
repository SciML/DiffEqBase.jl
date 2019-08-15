using DiffEqBase, Test
using Distributions

@test DiffEqBase.promote_tspan((0.0,1.0)) == (0.0,1.0)
@test DiffEqBase.promote_tspan((0,1.0)) == (0.0,1.0)
@test DiffEqBase.promote_tspan(1.0) == (0.0,1.0)
@test DiffEqBase.promote_tspan(nothing) == (nothing,nothing)

prob = ODEProblem((u,p,t)->u,(p,t0)->p[1],(p)->(0.0,p[2]),(2.0,1.0))
prob2 = DiffEqBase.get_concrete_problem(prob,nothing)

@test prob2.u0 == 2.0
@test prob2.tspan == (0.0,1.0)

prob = ODEProblem((u,p,t)->u,(p,t)->Normal(p,1),(0.0,1.0),1.0)
prob2 = DiffEqBase.get_concrete_problem(prob,nothing)
@test typeof(prob2.u0) == Float64

kwargs(; kw...) = kw
prob = ODEProblem((u,p,t)->u,1.0,nothing)
prob2 = DiffEqBase.get_concrete_problem(prob,kwargs(tspan=(1.2,3.4)))
@test prob2.tspan === (1.2,3.4)

prob = ODEProblem((u,p,t)->u,nothing,nothing)
prob2 = DiffEqBase.get_concrete_problem(prob,kwargs(u0=1.01,tspan=(1.2,3.4)))
@test prob2.u0 === 1.01

prob = ODEProblem((u,p,t)->u,1.0,(0,1))
@test_logs (:warn,
            "Integer time values are incompatible with adaptive integrators. Utilize " *
            "floating point numbers instead of integers in this case, i.e. (0.0,1.0) " *
            "instead of (0,1).") DiffEqBase.adaptive_warn(prob.u0,prob.tspan)

prob = DDEProblem((u, h, p, t) -> - h(p, t - p[1]), (p, t0) -> p[2], (p, t) -> 0,
                  (p) -> (0.0, p[3]), (1.0, 2.0, 3.0); constant_lags = (p) -> [p[1]])
prob2 = DiffEqBase.get_concrete_problem(prob, nothing)

@test prob2.u0 == 2.0
@test prob2.tspan == (0.0, 3.0)
@test prob2.constant_lags == [1.0]
