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

prob = ODEProblem((u,p,t)->u,1.0,(0,1))
DiffEqBase.adaptive_warn(prob.u0,prob.tspan)
