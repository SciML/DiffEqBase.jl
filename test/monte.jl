using StochasticDiffEq, DiffEqBase, DiffEqProblemLibrary, OrdinaryDiffEq
using Base.Test


prob = prob_sde_2Dlinear
prob2 = MonteCarloProblem(prob)
sim = solve(prob2,SRIW1(),dt=1//2^(3),num_monte=10)
calculate_monte_errors(sim)

prob = prob_sde_additivesystem
prob2 = MonteCarloProblem(prob)
sim = solve(prob2,SRA1(),dt=1//2^(3),num_monte=10)
calculate_monte_errors(sim)

output_func = function (x)
  last(last(x))^2
end
prob2 = MonteCarloProblem(prob,output_func=output_func)
sim = solve(prob2,SRA1(),dt=1//2^(3),num_monte=10)

prob = prob_sde_lorenz
prob2 = MonteCarloProblem(prob)
sim = solve(prob2,SRIW1(),dt=1//2^(3),num_monte=10)

prob = prob_ode_linear
prob_func = function (prob,i)
  prob.u0 = rand()*prob.u0
  prob
end
prob2 = MonteCarloProblem(prob,prob_func=prob_func)
sim = solve(prob2,Tsit5(),num_monte=100)
