using OrdinaryDiffEq, RecursiveArrayTools, Test, StaticArrays, DiffEqCallbacks, ForwardDiff


f = function (u,p,t)
  - u + sin(-t)
end


prob = ODEProblem(f,1.0,(0.0,-10.0))

condition= function (u,t,integrator) # Event when event_f(u,t,k) == 0
  - t - 2.95
end

affect! = function (integrator)
  integrator.u = integrator.u + 2
end

callback = ContinuousCallback(condition,affect!)

sol = solve(prob,Tsit5(),callback=callback)
@test length(sol) < 20

condition= function (out, u,t,integrator) # Event when event_f(u,t,k) == 0
  out[1] = - t - 2.95
end

affect! = function (integrator, idx)
  if idx == 1
    integrator.u = integrator.u + 2
  end
end

callback = VectorContinuousCallback(condition,affect!,1)

sol = solve(prob,Tsit5(),callback=callback)

f = function (du,u,p,t)
  du[1] = - u[1] + sin(t)
end


prob = ODEProblem(f,[1.0],(0.0,10.0))

condition= function (u,t,integrator) # Event when event_f(u,t,k) == 0
  t - 2.95
end

affect! = function (integrator)
  integrator.u = integrator.u .+ 2
end

callback = ContinuousCallback(condition,affect!)

sol = solve(prob,Tsit5(),callback=callback,abstol=1e-8,reltol=1e-6)

f = function (du,u,p,t)
  du[1] = u[2]
  du[2] = -9.81
end

condition= function (u,t,integrator) # Event when event_f(u,t,k) == 0
  u[1]
end

affect! = nothing
affect_neg! = function (integrator)
  integrator.u[2] = -integrator.u[2]
end

callback = ContinuousCallback(condition,affect!,affect_neg!,interp_points=100)

u0 = [50.0,0.0]
tspan = (0.0,15.0)
prob = ODEProblem(f,u0,tspan)


sol = solve(prob,Tsit5(),callback=callback,adaptive=false,dt=1/4)

condition_single = function (u,t,integrator) # Event when event_f(u,t,k) == 0
  u
end

affect! = nothing
affect_neg! = function (integrator)
  integrator.u[2] = -integrator.u[2]
end

callback_single = ContinuousCallback(condition_single,affect!,affect_neg!,interp_points=100,idxs=1)

u0 = [50.0,0.0]
tspan = (0.0,15.0)
prob = ODEProblem(f,u0,tspan)

sol = solve(prob,Tsit5(),callback=callback_single,adaptive=false,dt=1/4)
sol = solve(prob,Tsit5(),callback=callback_single,save_everystep=false)
t = sol.t[end÷2] # this is the callback time point
sol = solve(prob,Tsit5(),callback=callback_single,saveat=t)
@test count(x->x==t, sol.t) == 2
sol = solve(prob,Tsit5(),callback=callback_single,saveat=t-eps(t))
@test count(x->x==t, sol.t) == 2
# check interpolation @ discontinuity
@test sol(t,continuity=:right)[2] > 0
@test sol(t,continuity=:left)[2] < 0

#plot(sol,denseplot=true)

sol = solve(prob,Vern6(),callback=callback)
#plot(sol,denseplot=true)
sol = solve(prob,BS3(),callback=callback)

sol33 = solve(prob,Vern7(),callback=callback)

bounced = ODEProblem(f,sol[8],(0.0,1.0))
sol_bounced = solve(bounced,Vern6(),callback=callback,dt=sol.t[9]-sol.t[8])
#plot(sol_bounced,denseplot=true)
sol_bounced(0.04) # Complete density
@test maximum(maximum.(map((i)->sol.k[9][i]-sol_bounced.k[2][i],1:length(sol.k[9])))) == 0


sol2= solve(prob,Vern6(),callback=callback,adaptive=false,dt=1/2^4)
#plot(sol2)

sol2= solve(prob,Vern6())

sol3= solve(prob,Vern6(),saveat=[.5])

## Saving callback

condition = function (u,t,integrator)
  true
end
affect! = function (integrator) end

save_positions = (true,false)
saving_callback = DiscreteCallback(condition,affect!,save_positions=save_positions)

sol4 = solve(prob,Tsit5(),callback=saving_callback)

@test sol2(3) ≈ sol(3)

affect! = function (integrator)
  u_modified!(integrator,false)
end
saving_callback2 = DiscreteCallback(condition,affect!,save_positions=save_positions)
sol4 = solve(prob,Tsit5(),callback=saving_callback2)

cbs = CallbackSet(saving_callback,saving_callback2)
sol4_extra = solve(prob,Tsit5(),callback=cbs)

@test length(sol4_extra) == 2length(sol4) - 1

condition= function (u,t,integrator)
  u[1]
end

affect! = function (integrator, retcode = nothing)
  if retcode === nothing
    terminate!(integrator)
  else
    terminate!(integrator, retcode)
  end
end

terminate_callback = ContinuousCallback(condition,affect!)
custom_retcode_callback = ContinuousCallback(condition,x->affect!(x,:Custom))

tspan2 = (0.0,Inf)
prob2 = ODEProblem(f,u0,tspan2)

sol5 = solve(prob2,Tsit5(),callback=terminate_callback)
sol5_1 = solve(prob2,Tsit5(),callback=custom_retcode_callback)

@test sol5.retcode == :Terminated
@test sol5_1.retcode == :Custom
@test sol5[end][1] < 3e-12
@test sol5.t[end] ≈ sqrt(50*2/9.81)

affect2! = function (integrator)
  if integrator.t >= 3.5
    terminate!(integrator)
  else
    integrator.u[2] = -integrator.u[2]
  end
end
terminate_callback2 = ContinuousCallback(condition,nothing,affect2!,interp_points=100)


sol5 = solve(prob2,Vern7(),callback=terminate_callback2)

@test sol5[end][1] < 1.3e-10
@test sol5.t[end] ≈ 3*sqrt(50*2/9.81)

condition= function (u,t,integrator) # Event when event_f(u,t,k) == 0
  t-4
end

affect! = function (integrator)
  terminate!(integrator)
end

terminate_callback3 = ContinuousCallback(condition,affect!,interp_points=1000)

bounce_then_exit = CallbackSet(callback,terminate_callback3)

sol6 = solve(prob2,Vern7(),callback=bounce_then_exit)

@test sol6[end][1] > 0
@test sol6[end][1] < 100
@test sol6.t[end] ≈ 4

# Test ContinuousCallback hits values on the steps
t_event = 100.0
f_simple(u,p,t) = 1.00001*u
event_triggered = false
condition_simple(u,t,integrator) = t_event-t
function affect_simple!(integrator)
  global event_triggered
  event_triggered = true
end
cb = ContinuousCallback(condition_simple, nothing, affect_simple!)
prob = ODEProblem(f_simple, [1.0], (0.0, 2.0*t_event))
sol = solve(prob,Tsit5(),callback=cb, adaptive = false, dt = 10.0)
@test event_triggered

# https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl/issues/328
ode = ODEProblem((du, u, p, t) -> (@. du .= -u), ones(5), (0.0, 100.0))
sol = solve(ode, AutoTsit5(Rosenbrock23()), callback=TerminateSteadyState())
sol1 = solve(ode, Tsit5(), callback=TerminateSteadyState())
@test sol.u == sol1.u

# DiscreteCallback
f = function (du,u,p,t)
  du[1] = -0.5*u[1] + 10
  du[2] = -0.5*u[2]
end

u0 = [10,10.]
tstop = [5.;8.]
prob = ODEProblem(f, u0, (0, 10.))
condition = (u,t,integrator) -> t in tstop
affect! = (integrator) -> integrator.u .= 1.0
save_positions = (true,true)
cb = DiscreteCallback(condition, affect!, save_positions=save_positions)
sol1 = solve(prob,Tsit5(),callback = cb,tstops=tstop,saveat=tstop)
@test count(x->x==tstop[1], sol1.t) == 2
@test count(x->x==tstop[2], sol1.t) == 2
sol2 = solve(prob,Tsit5(),callback = cb,tstops=tstop,saveat=prevfloat.(tstop))
@test count(x->x==tstop[1], sol2.t) == 2
@test count(x->x==tstop[2], sol2.t) == 2

function model(du, u, p, t)
    du[1] = 0.
    for i in 2:(length(du)-1)
        du[i] = p[i] * (u[i-1] - u[i])
    end
    du[end] = p[end] * (p[1] * u[end-1] - u[end])
    return nothing
end

perror = [1.0, 0.02222434508140991, 0.017030281542289794, 0.015917011145559996, 0.1608874463597176, 0.13128016561792297, 0.11056834258380167, 0.5222141958458832, 1.0711942201995688, 0.2672878398678257, 8.900058706990183, 0.010760065201065117, 0.016319181296867765, 2.2693845639611925, 0.2152216345154439, 0.029186712540925457, 0.21419429135100806, 0.029177617589788596, 0.03064986043089549, 0.023280222517122397, 6.931251277770224]
y_max = 0.002604806609572015
u0 = [1, zeros(length(perror) - 1)...]
tspan = (0., 5000.)

condition = (u, t, i) -> (t == 1.)
affect! = (i) -> (i.u[1] = 0.)

condition2 = (u, t, i) -> (u[end] - y_max / 2.)
t_half_1 = 0.
affect2! = function (i)
  global t_half_1 = i.t
end

prob = ODEProblem(model, u0, tspan, perror)
integrator = init(
    prob,
    Rosenbrock23();
    callback=CallbackSet(
        PositiveDomain(),
        DiscreteCallback(condition, affect!),
        ContinuousCallback(condition2, affect2!, terminate!),
    ),
    tstops = [1.],
    force_dtmin=true,
    progress=true
)

sol = solve!(integrator)


### https://github.com/SciML/DifferentialEquations.jl/issues/662

function prob1!(du, u, p, s)
  du[1] = u[2]
  du[2] = zero(eltype(u))
end

function cb1test(out, u, x)
  out[1] = u[1] - x
end

function cb1action!(i)
  i.u[2] = -i.u[2]
end

function ode1(x)
  prob = ODEProblem{true}(prob1!, [zero(typeof(x)), one(typeof(x))], (0.0, 1.0))
  cb = VectorContinuousCallback(
    (out, u, t, i) -> cb1test(out, u, x),
    (i, ndx) -> cb1action!(i),
    1; rootfind=true)
  soln = solve(prob, Tsit5(); callback=cb)
  soln[end][1]
end

@test ForwardDiff.derivative(ode1, 0.42) ≈ 2
