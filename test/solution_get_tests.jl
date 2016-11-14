using OrdinaryDiffEq, StochasticDiffEq, DiffEqDevTools, DiffEqProblemLibrary
prob = prob_ode_2Dlinear

## Solve and plot
println("Test getindex")
sol =solve(prob,abstol=1e-8,reltol=1e-7)

length(sol)
sol[1]
sol[1,2]
print(sol)
show(sol)

srand(100)
dts = 1./2.^(10:-1:4) #14->7 good plot

prob2 = prob_sde_wave
sim = test_convergence(dts,prob2,numMonte=Int(1e1),alg=:EM)

length(sim)
sim[1]
sim[1][1]
print(sim)
show(sim)
sim[end]

sim = monteCarloSim(prob2,dt=1/2^(3),save_timeseries=true,numMonte=5)
length(sim)
sim[1]
sim[1][1]
print(sim)
show(sim)
sim[end]

@test sol[1] == prob.u0
