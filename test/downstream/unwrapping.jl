using OrdinaryDiffEq, Test
my_f(u, p, t) = u
my_f!(du, u, p, t) = du .= u
ode = ODEProblem(my_f, [1.0], (0.0, 1.0))
integrator = init(ode, Tsit5())
@test SciMLBase.unwrapped_f(integrator.f.f) === my_f

ode = ODEProblem(my_f!, [1.0], (0.0, 1.0))
integrator = init(ode, Tsit5())
@test SciMLBase.unwrapped_f(integrator.f.f) === my_f!