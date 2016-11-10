using OrdinaryDiffEq, DiffEqProblemLibrary, DiffEqBase, Sundials

alg, kwargs = default_algorithm(prob_ode_2Dlinear;dt=1//2^(4))
sol =solve(prob_ode_2Dlinear;dt=1//2^(4))

bool1 = typeof(sol.alg) == alg
bool2 = typeof(sol.alg) == Tsit5

sol =solve(prob_ode_2Dlinear;alg_hints=[:stiff])

bool3 = sol.alg == CVODE_BDF

const linear_bigα = parse(BigFloat,"1.01")
f = (t,u,du) -> begin
  for i in 1:length(u)
    du[i] = linear_bigα*u[i]
  end
end
analytic = (t,u0) -> u0*exp(linear_bigα*t)
prob_ode_bigfloat2Dlinear = ODETestProblem(f,map(BigFloat,rand(4,2)).*ones(4,2)/2,analytic)

sol =solve(prob_ode_bigfloat2Dlinear;dt=1//2^(4))
bool4 = typeof(sol.alg) == Vern8

default_algorithm(prob_ode_bigfloat2Dlinear;alg_hints=[:stiff])

sol =solve(prob_ode_bigfloat2Dlinear;alg_hints=[:stiff])

bool5 = typeof(sol.alg) == Rosenbrock23

bool1 && bool2 && bool3 && bool4 && bool5
