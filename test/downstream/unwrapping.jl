using OrdinaryDiffEq, Test
my_f(u, p, t) = u
my_f!(du, u, p, t) = du .= u
ode = ODEProblem(my_f, [1.0], (0.0, 1.0))
integrator = init(ode, Tsit5())
@test SciMLBase.unwrapped_f(integrator.f.f) === my_f

ode = ODEProblem(my_f!, [1.0], (0.0, 1.0))
integrator = init(ode, Tsit5())
@test SciMLBase.unwrapped_f(integrator.f.f) === my_f!

using ForwardDiff, Measurements
x = 1.0 ± 0.0
f = (du,u,p,t)-> du .= u
tspan = (ForwardDiff.Dual(0.0,(0.01)),ForwardDiff.Dual(1.0,(0.01)))
prob = ODEProblem(f, [x], tspan)

# Should not error during problem construction but should be unwrapped
integ = init(prob, Tsit5(), dt = 0.1)
@test integ.f.f === f