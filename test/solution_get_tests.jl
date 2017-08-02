using OrdinaryDiffEq, StochasticDiffEq, DiffEqDevTools, DiffEqProblemLibrary
prob = prob_ode_2Dlinear

## Solve and plot
println("Test getindex")
sol = solve(prob,Tsit5(),abstol=1e-8,reltol=1e-7)

length(sol)
sol[1]
sol[1,2]
#print(sol)
#show(sol)

srand(100)
dts = 1./2.^(10:-1:4) #14->7 good plot

prob2 = prob_sde_wave
sol = test_convergence(dts,prob2,EM(),numMonte=Int(1e1))

length(sol)
sol[1]
sol[1][1]
sol[end]

@test sol[1] == prob2.u0
