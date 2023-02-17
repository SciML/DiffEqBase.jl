using OrdinaryDiffEq

# Auto callback merging

do_nothing = DiscreteCallback((u, t, integrator) -> true,
                              integrator -> nothing)
problem = ODEProblem((u, p, t) -> -u,
                     1.0, (0.0, 1.0),
                     callback = do_nothing)
solve(problem, Euler(),
      dt = 0.1,
      callback = do_nothing)
