using OrdinaryDiffEq, ForwardDiff, Zygote, Test
using SciMLSensitivity
using Random

function dt!(du, u, p, t)
    x, y = u
    α, β, δ, γ = p
    du[1] = dx = α * x - β * x * y
    du[2] = dy = -δ * y + γ * x * y
end

n_par = 3
Random.seed!(2)
u0 = rand(2, n_par)
u0[:, 1] = [1.0, 1.0]
tspan = (0.0, 10.0)
p = [2.2, 1.0, 2.0, 0.4]
prob_ode = ODEProblem(dt!, u0[:, 1], tspan)

function test_loss(p1, prob)
    function prob_func(prob, i, repeat)
        @show i
        remake(prob, u0 = u0[:, i])
    end

    #define ensemble problem
    ensembleprob = EnsembleProblem(prob, prob_func = prob_func)

    u = Array(solve(ensembleprob, Tsit5(), EnsembleThreads(), trajectories = n_par,
                    p = p1,
                    sensealg = ForwardDiffSensitivity(),
                    saveat = 0.1, dt = 0.001))[:, end, :]
    loss = sum(u)
    return loss
end

test_loss(p, prob_ode)

@time gs = Zygote.gradient(p) do p
    test_loss(p, prob_ode)
end
@test gs[1] isa Vector

### https://github.com/SciML/DiffEqFlux.jl/issues/595

function fiip(du, u, p, t)
    du[1] = dx = p[1] * u[1] - p[2] * u[1] * u[2]
    du[2] = dy = -p[3] * u[2] + p[4] * u[1] * u[2]
end

p = [1.5, 1.0, 3.0, 1.0];
u0 = [1.0; 1.0];
prob = ODEProblem(fiip, u0, (0.0, 10.0), p)
sol = solve(prob, Tsit5())

function sum_of_solution(x)
    _prob = remake(prob, u0 = x[1:2], p = x[3:end])
    sum(solve(_prob, Tsit5(), saveat = 0.1))
end
Zygote.gradient(sum_of_solution, [u0; p])

# Testing ensemble problem. Works with ForwardDiff. Does not work with Zygote.
N = 3
eu0 = rand(N, 2)
ep = rand(N, 4)

ensemble_prob = EnsembleProblem(prob,
                                prob_func = (prob, i, repeat) -> remake(prob,
                                                                        u0 = eu0[i, :],
                                                                        p = ep[i, :],
                                                                        saveat = 0.1))
esol = solve(ensemble_prob, Tsit5(), trajectories = N)

cache = Ref{Any}()

function sum_of_e_solution(p)
    ensemble_prob = EnsembleProblem(prob,
                                    prob_func = (prob, i, repeat) -> remake(prob,
                                                                            u0 = eu0[i, :],
                                                                            p = p[i, :],
                                                                            saveat = 0.1))
    sol = solve(ensemble_prob, Tsit5(), EnsembleSerial(), trajectories = N, abstol = 1e-12,
                reltol = 1e-12)
    z = Array(sol[1])
    cache[] = sol[1].t
    sum(z) # just test for the first solutions, gradients should be zero for others
end

sum_of_e_solution(ep)

x = ForwardDiff.gradient(sum_of_e_solution, ep)
y = Zygote.gradient(sum_of_e_solution, ep)[1] # Zygote second to test cache of forward pass
@test x ≈ y
@test cache[] == 0:0.1:10.0 # test prob.kwargs is forwarded
