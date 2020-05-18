using OrdinaryDiffEq, Test
function lorenz(du,u,p,t)
 du[1] = 10.0(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - (8/3)*u[3]
end
u0 = [1.0;0.0;0.0]
tspan = (0.0,1.0)
prob = ODEProblem(lorenz,u0,tspan)
sol = solve(prob,Tsit5())
@test_broken @inferred solve(prob,Tsit5())

function f(du,u,p,t)
  du[1] = p.a
  du[2] = p.b
end

const alg = Tsit5()

function solve_ode(f::F, p::P) where {F,P}

  tspan = (0., 1.0)
  Δt = tspan[2] - tspan[1]
  dt = 1/252
  nodes = Int(ceil(Δt / dt) + 1)
  t = T = [tspan[1] + (i - 1) * dt for i = 1:nodes]

  # if I do not set {true}, prob type Any...
  prob = ODEProblem{true}(f, [0., 0.], tspan, p)
  # prob = ODEProblem(f, [0., 0.], tspan, p)

  prob_func = (prob, i, repeat) -> begin
    remake(prob, tspan = (T[i + 1], t[1]))
  end

  # ensemble problem
  odes = EnsembleProblem(prob, prob_func = prob_func)

  sol = OrdinaryDiffEq.solve(
    odes, OrdinaryDiffEq.Tsit5(), OrdinaryDiffEq.EnsembleThreads(),
    trajectories = nodes - 1, saveat = -dt
  )

  return sol
end
@test_broken @inferred solve_ode(f, (a = 1, b = 1))

using StochasticDiffEq, Test
u0=1/2
ff(u,p,t) = u
gg(u,p,t) = u
dt = 1//2^(4)
tspan = (0.0,1.0)
prob = SDEProblem(ff,gg,u0,(0.0,1.0))
sol = solve(prob,EM(),dt=dt)
@test_broken @inferred solve(prob,EM(),dt=dt)
