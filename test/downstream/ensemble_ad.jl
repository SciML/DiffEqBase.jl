using OrdinaryDiffEq, Zygote
using DiffEqSensitivity
using Random

function dt!(du, u, p, t)
  x, y = u
  α, β, δ, γ = p
  du[1] = dx = α*x - β*x*y
  du[2] = dy = -δ*y + γ*x*y
end

n_par=3
Random.seed!(2)
u0=rand(2,n_par)
u0[:,1] = [1.0,1.0]
tspan = (0.0, 10.0)
p = [2.2, 1.0, 2.0, 0.4]
prob_ode = ODEProblem(dt!, u0[:,1], tspan)

function test_loss(p1,prob)

	function prob_func(prob, i, repeat)
		@show i
		remake(prob,u0=u0[:,i])
	end

	#define ensemble problem
	ensembleprob = EnsembleProblem(prob,prob_func = prob_func)

	u = Array(solve(ensembleprob, Tsit5(), EnsembleThreads(), trajectories=n_par,
	p=p,
	sensealg = ForwardDiffSensitivity(),
	saveat = 0.1, dt=0.001))[:,end,:]
	loss=sum(u)
	return loss
end

@time gs = Zygote.gradient(p) do p
	test_loss(p,prob_ode)
end
