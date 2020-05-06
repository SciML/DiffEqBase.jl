using OrdinaryDiffEq, StochasticDiffEq, SteadyStateDiffEq, RecursiveArrayTools, Test

function f!(du,u,p,t)
  du[1] = p[1] + p[2]*u[1]
  du[2] = p[3]*u[1] + p[4]*u[2]
end
u0 = zeros(2)
p = [2.0,-2.0,1.0,-4.0]

probODE = ODEProblem(f!,u0,(0.0,10.0),p)
probSS = SteadyStateProblem(f!,u0,p)

solODE = concrete_solve(probODE,Tsit5(), abstol=1e-14,reltol=1e-14)
solSS = concrete_solve(probSS,DynamicSS(Rodas5()))

@test solSS ≈ solODE.u[end] rtol = 1e-8
@test length(fieldnames(typeof(solODE))) == 2
@test length(fieldnames(typeof(solSS))) == 0
