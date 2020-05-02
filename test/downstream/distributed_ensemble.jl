using Distributed
addprocs(2)
println("There are $(nprocs()) processes")
@everywhere using OrdinaryDiffEq

@everywhere prob = ODEProblem((u,p,t)->1.01u,0.5,(0.0,1.0))
@everywhere u0s = [rand()*prob.u0 for i in 1:2]
@everywhere function prob_func(prob,i,repeat)
    println("Running trajectory $i")
    ODEProblem(prob.f,u0s[i],prob.tspan)
end

ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)
sim = solve(ensemble_prob,Tsit5(),EnsembleSplitThreads(),trajectories=2)
