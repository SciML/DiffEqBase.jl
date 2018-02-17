using DiffEqBase: @add_kwonly, add_kwonly

@add_kwonly function f(a, b; c=3, d=4)
  (a, b, c, d)
end
@test f(1, 2) == (1, 2, 3, 4)
@test f(a=1, b=2) == (1, 2, 3, 4)
@test_throws ErrorException f()

@add_kwonly g(a, b; c=3, d=4) = (a, b, c, d)
@test g(1, 2) == (1, 2, 3, 4)
@test g(a=1, b=2) == (1, 2, 3, 4)

@add_kwonly h(; c=3, d=4) = (c, d)
@test h() == (3, 4)

@test_throws ErrorException add_kwonly(:(i(c=3, d=4) = (c, d)))


function f(du,u,p,t)
  du[1] = 0.2u[1]
  du[2] = 0.4u[2]
end
u0 = ones(2)
tspan = (0,1.0)

# Create a ODEProblem and test set_u0:
prob1 = SplitODEProblem(f,f,u0,tspan,Dict(); mass_matrix=eye(length(u0)))
prob2 = remake(prob1; u0 = prob1.u0 + 1)
@test prob1.f === prob2.f
@test prob1.p === prob2.p
@test prob1.u0 + 1 ≈ prob2.u0
@test prob1.tspan == prob2.tspan
@test prob1.jac_prototype === prob2.jac_prototype
@test prob1.callback === prob2.callback
@test prob1.mass_matrix === prob2.mass_matrix
@test prob1.problem_type === prob2.problem_type
