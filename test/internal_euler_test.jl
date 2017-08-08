using DiffEqBase, DiffEqBase.InternalEuler

# Try it
ff(t,u) = u
ff(::Val{:jac}, t, u) = 1
ff(::Type{Val{:analytic}},t,u0) = u0*exp(t)

dt = 0.01
prob = ODEProblem(ff, 1.0, (0.0,1.0) )
sol = solve(prob, InternalEuler.FwdEulerAlg(), tstops=0:dt:1)
sol2 = solve(prob, InternalEuler.BwdEulerAlg(), tstops=0:dt:1)

#using Plots
#plot(sol)
#plot!(sol2)
#plot!(sol2, plot_analytic=true)

#using DiffEqDevTools
#dts = 1./2.^(8:-1:4)
#sim = test_convergence(dts,p2,BwdEulerAlg())
#@show sim.𝒪est[:final]
#plot(sim)
