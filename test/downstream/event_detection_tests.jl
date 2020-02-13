using StaticArrays
using OrdinaryDiffEq
using Test

@inbounds @inline function ż(z, p, t)
    A, B, D = p
    p₀, p₂ = z[1:2]
    q₀, q₂ = z[3:4]

    return SVector{4}(
        -A * q₀ - 3 * B / √2 * (q₂^2 - q₀^2) - D * q₀ * (q₀^2 + q₂^2),
        -q₂ * (A + 3 * √2 * B * q₀ + D * (q₀^2 + q₂^2)),
        A * p₀,
        A * p₂
  )
end

condition(u, t, integrator) = u
affect!(integrator) = nothing
cbf(idx) = ContinuousCallback(condition,
    affect!, nothing, save_positions=(false, true), idxs=idx)
z0 = SVector{4}(7.1989885061904335, -0.165912283356219, 0., -3.63534900748947)

tspan=(0.,300.)
prob=ODEProblem(ż, z0, tspan, (A=1, B=0.55, D=0.4), callback=cbf(3))
sol=solve(prob, Vern9(), abstol=1e-14, reltol=1e-14,
    save_everystep=false, save_start=false, save_end=false, maxiters=1e6)

@test length(sol) > 100

prob=ODEProblem(ż, z0, (0,400.), (A=1, B=0.55, D=0.4), callback=cbf(3))
sol=solve(prob, Vern9(), abstol=1e-14, reltol=1e-14, save_everystep=false, save_start=false, save_end=false, maxiters=2e4)

@test length(sol) > 100

prob=ODEProblem(ż, z0, (0,5000.), (A=1, B=0.55, D=0.4), callback=cbf(3))
sol=solve(prob, Vern9(), abstol=1e-14, reltol=1e-14, save_everystep=false, save_start=false, save_end=false, maxiters=1e6)

@test length(sol) > 1500

f = function (du,u,p,t)
  du[1] = u[2]
  du[2] = -p[1]
end
function condition(u,t,integrator) # Event when event_f(u,t) == 0
  u[1]
end
function affect!(integrator)
  @test integrator.u[1] >= 0
  integrator.u[2] = -integrator.u[2]
end
cb2 = ContinuousCallback(condition,affect!)
tspan = (0.0,10000.0)
u0 = [50.0,0.0]
p = 9.8
prob = ODEProblem(f,u0,tspan,p)
sol = solve(prob,Tsit5(),callback=cb2)
@test minimum(sol') > -40

function vcondition!(out,u,t,integrator)
  out[1] = u[1]
  out[2] = u[2]
end

function vaffect!(integrator, event_idx)
  @test integrator.u[1] >= 0.0
  if event_idx == 1
    integrator.u[2] = -integrator.u[2]
  else
    integrator.p = 0.0
  end
end

u0 = [50.0,0.0]
tspan = (0.0,15.0)
p = 9.8
prob = ODEProblem(f,u0,tspan,p)
Vcb = VectorContinuousCallback(vcondition!,vaffect!, 2 , save_positions=(true,true))
sol = solve(prob,Tsit5(), callback=Vcb)

f = function (du,u,p,t)
  du[1] = u[2]
  du[2] = -9.81
end

function condition(u,t,integrator) # Event when event_f(u,t) == 0
  u[1]    # Event when height crosses from positive to negative
end

function affect_neg(integrator)
  @test integrator.u[1] >= 0.0
  integrator.u[2] = -0.8*integrator.u[2]
end

cb = ContinuousCallback(condition,nothing,affect_neg! = affect_neg)

u0 = [1.0,0.0]
tspan = (0.0, 3.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob,Tsit5(),saveat=0.01,callback=cb)

f! = function (du,u,p,t)
  du[1] = u[2]
  du[2] = -p
end

function condition!(out,u,t,integrator) 
  out[1] = u[1]
  out[2] = u[2]
end

function affect!(integrator, event_idx)
  if event_idx == 1
    integrator.u[2] = -integrator.u[2]
  else
    integrator.p = 0.0
  end
end

u0 = [50.0,0.0]
tspan = (0.0,15.0)

begin
    p = 9.8
    prob = ODEProblem(f!,u0,tspan,p)
    Vcb = VectorContinuousCallback(condition!,affect!, 2 , save_positions=(true,true))
    sol = solve(prob,Tsit5(), callback=Vcb)
end

