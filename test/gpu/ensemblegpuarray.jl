using DiffEqGPU, CuArrays, OrdinaryDiffEq, Test

function lorenz(du,u,p,t)
@inbounds begin
   du[1] = p[1]*(u[2]-u[1])
   du[2] = u[1]*(p[2]-u[3]) - u[2]
   du[3] = u[1]*u[2] - p[3]*u[3]
end
nothing
end

CuArrays.allowscalar(false)
u0 = Float32[1.0;0.0;0.0]
tspan = (0.0f0,100.0f0)
p = (10.0f0,28.0f0,8/3f0)
prob = ODEProblem(lorenz,u0,tspan,p)
const pre_p = [rand(Float32,3) for i in 1:10]
prob_func = (prob,i,repeat) -> remake(prob,p=pre_p[i].*p)
monteprob = EnsembleProblem(prob, prob_func = prob_func)

@info "Explicit Methods"

#Performance check with nvvp
# CUDAnative.CUDAdrv.@profile
@time sol = solve(monteprob,Tsit5(),EnsembleGPUArray(),trajectories=10,saveat=1.0f0)
@test length(filter(x -> x.u != sol.u[1].u, sol.u)) != 0 # 0 element array
@time sol = solve(monteprob,ROCK4(),EnsembleGPUArray(),trajectories=10,saveat=1.0f0)
@time sol2 = solve(monteprob,Tsit5(),EnsembleGPUArray(),trajectories=10,
                                               batch_size=5,saveat=1.0f0)

@test length(filter(x -> x.u != sol.u[1].u, sol.u)) != 0 # 0 element array
@test length(filter(x -> x.u != sol2.u[6].u, sol.u)) != 0 # 0 element array
@test all(all(sol[i].prob.p .== pre_p[i].*p) for i in 1:10)
@test all(all(sol2[i].prob.p .== pre_p[i].*p) for i in 1:10)

@time solve(monteprob,Tsit5(),EnsembleCPUArray(),trajectories=10,saveat=1.0f0)
@time solve(monteprob,Tsit5(),EnsembleThreads(), trajectories=10,saveat=1.0f0)
@time solve(monteprob,Tsit5(),EnsembleSerial(),  trajectories=10,saveat=1.0f0)

#=
solve(monteprob,TRBDF2(),EnsembleCPUArray(),dt=0.1,trajectories=2,saveat=1.0f0)
solve(monteprob,TRBDF2(),EnsembleGPUArray(),dt=0.1,trajectories=2,saveat=1.0f0)
@test_broken solve(monteprob,TRBDF2(linsolve=LinSolveGPUSplitFactorize()),EnsembleGPUArray(),dt=0.1,trajectories=2,saveat=1.0f0)
=#

@info "Implicit Methods"

function lorenz_jac(J,u,p,t)
@inbounds begin
   σ = p[1]
   ρ = p[2]
   β = p[3]
   x = u[1]
   y = u[2]
   z = u[3]
   J[1,1] = -σ
   J[2,1] = ρ - z
   J[3,1] = y
   J[1,2] = σ
   J[2,2] = -1
   J[3,2] = x
   J[1,3] = 0
   J[2,3] = -x
   J[3,3] = -β
end
nothing
end

function lorenz_tgrad(J,u,p,t)
nothing
end

func = ODEFunction(lorenz,jac=lorenz_jac,tgrad=lorenz_tgrad)
prob_jac = ODEProblem(func,u0,tspan,p)
monteprob_jac = EnsembleProblem(prob_jac, prob_func = prob_func)

@time solve(monteprob_jac,Rodas5(linsolve=LinSolveGPUSplitFactorize()),EnsembleCPUArray(),dt=0.1,trajectories=10,saveat=1.0f0)
@time solve(monteprob_jac,Rodas5(linsolve=LinSolveGPUSplitFactorize()),EnsembleGPUArray(),dt=0.1,trajectories=10,saveat=1.0f0)
@time solve(monteprob_jac,TRBDF2(linsolve=LinSolveGPUSplitFactorize()),EnsembleCPUArray(),dt=0.1,trajectories=10,saveat=1.0f0)
@time solve(monteprob_jac,TRBDF2(linsolve=LinSolveGPUSplitFactorize()),EnsembleGPUArray(),dt=0.1,trajectories=10,saveat=1.0f0)

@info "Callbacks"

condition = function (u,t,integrator)
  @inbounds u[1] > 5
end

affect! = function (integrator)
  @inbounds integrator.u[1] = -4
end

callback_prob = ODEProblem(lorenz,u0,tspan,p,callback=DiscreteCallback(condition,affect!,save_positions=(false,false)))
callback_monteprob = EnsembleProblem(callback_prob, prob_func = prob_func)
@time solve(callback_monteprob,Tsit5(),EnsembleGPUArray(),trajectories=10,saveat=1.0f0)

c_condition = function (u,t,integrator)
  @inbounds u[1] - 3
end

c_affect! = function (integrator)
  @inbounds integrator.u[1] += 20
end

callback_prob = ODEProblem(lorenz,u0,tspan,p,callback=ContinuousCallback(c_condition,c_affect!,save_positions=(false,false)))
callback_monteprob = EnsembleProblem(callback_prob, prob_func = prob_func)
solve(callback_monteprob,Tsit5(),EnsembleGPUArray(),trajectories=2,saveat=1.0f0)

@info "ROBER"

#=
using OrdinaryDiffEq, LinearAlgebra, ParameterizedFunctions
LinearAlgebra.BLAS.set_num_threads(1)
rober = @ode_def begin
dy₁ = -k₁*y₁+k₃*y₂*y₃
dy₂ =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
dy₃ =  k₂*y₂^2
end k₁ k₂ k₃
=#

function rober_f(internal_var___du,internal_var___u,internal_var___p,t)
  @inbounds begin
      internal_var___du[1] = -(internal_var___p[1]) * internal_var___u[1] + internal_var___p[3] * internal_var___u[2] * internal_var___u[3]
      internal_var___du[2] = (internal_var___p[1] * internal_var___u[1] - internal_var___p[2] * internal_var___u[2] ^ 2) - internal_var___p[3] * internal_var___u[2] * internal_var___u[3]
      internal_var___du[3] = internal_var___p[2] * internal_var___u[2] ^ 2
  end
  nothing
end

function rober_jac(internal_var___J,internal_var___u,internal_var___p,t)
  @inbounds begin
      internal_var___J[1, 1] = -(internal_var___p[1])
      internal_var___J[1, 2] = internal_var___p[3] * internal_var___u[3]
      internal_var___J[1, 3] = internal_var___p[3] * internal_var___u[2]
      internal_var___J[2, 1] = internal_var___p[1] * 1
      internal_var___J[2, 2] = -2 * internal_var___p[2] * internal_var___u[2] - internal_var___p[3] * internal_var___u[3]
      internal_var___J[2, 3] = -(internal_var___p[3]) * internal_var___u[2]
      internal_var___J[3, 1] = 0 * 1
      internal_var___J[3, 2] = 2 * internal_var___p[2] * internal_var___u[2]
      internal_var___J[3, 3] = 0 * 1
  end
  nothing
end

function rober_tgrad(J,u,p,t)
nothing
end

rober_prob = ODEProblem(ODEFunction(rober_f,jac=rober_jac,tgrad=rober_tgrad),
                      Float32[1.0,0.0,0.0],(0.0f0,1f5),(0.04f0,3f7,1f4))
sol = solve(rober_prob,Rodas5(),abstol=1f-8,reltol=1f-8)
sol = solve(rober_prob,TRBDF2(),abstol=1f-4,reltol=1f-1)
rober_monteprob = EnsembleProblem(rober_prob, prob_func = prob_func)

@time sol = solve(rober_monteprob,Rodas5(linsolve=LinSolveGPUSplitFactorize()),
                EnsembleGPUArray(),trajectories=10,saveat=1.0f0,abstol=1f-8,
                reltol=1f-8)
@time sol = solve(rober_monteprob,TRBDF2(linsolve=LinSolveGPUSplitFactorize()),
                EnsembleGPUArray(),trajectories=10,saveat=1.0f0,abstol=1f-4,
                reltol=1f-1)
@time sol = solve(rober_monteprob,TRBDF2(),EnsembleThreads(),trajectories=10,
                abstol=1e-4,reltol=1e-1,saveat=1.0f0)
