using StaticArrays
using OrdinaryDiffEq
using LinearAlgebra
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


# https://github.com/SciML/DifferentialEquations.jl/issues/758
"""
    eoms(u, p, t)

Equations of motion for the system (6-Dimensional).
"""
function eoms(u, p, t)
    mu = p
    omm = 1 - mu

    r13 = sqrt((u[1] + mu)^2 + u[2]^2 + u[3]^2)
    r13_3 = r13^3
    r23 = sqrt((u[1] - 1 + mu)^2 + u[2]^2 + u[3]^2)
    r23_3 = r23^3

    Ωx = u[1] - (omm) * (u[1] + mu) / r13_3 - mu * (u[1] - omm) / r23_3
    Ωy = u[2] - (omm) * u[2] / r13_3 - mu * u[2] / r23_3
    Ωz = -(omm) * u[3] / r13_3 - mu * u[3] / r23_3

    return @SVector [
        u[4],
        u[5],
        u[6],
        2u[5] + Ωx,
        -2u[4] + Ωy,
        Ωz
    ]
end

"""
    distance(sol, q, t)

Calculate the isochronous position variation between `q` and `sol(t)`.
"""
function distance(sol, q, t)
    posinds = @SVector [1, 2, 3]
    q_sol = sol(t)
    return norm(q[posinds] .- q_sol[posinds])
end

"""
    distance_callback(sol, threshold)

Create a `ContinuousCallback` that terminates integration when `distance` is
less than `threshold`.
"""
function distance_callback(sol, threshold)
    condition(u, t, _) = distance(sol, u, t) - threshold
    affect!(integrator) = nothing
    affect_neg!(integrator) = terminate!(integrator)
    return ContinuousCallback(condition, affect!, affect_neg!)
end

# Mu Parameter
μ = 0.012150607114626023

# State of Chief at time 0
q_chief = @SVector [1.0220277279709828, 0.0, -0.18210117430516676, 0.0, -0.1032699633021813, 0.0]
tof_chief = 10.0

# Integrate the chief trajectory
prob_chief = ODEProblem{false}(eoms, q_chief, (0.0, tof_chief), μ)
sol_chief = solve(prob_chief, Vern9(); abstol=1e-12, reltol=1e-12)

# Define point where deputy departs the chief
t_dep_start = 0.012113772215335508 # Time the deputy trajectory starts at

# Get state of chief at this departure point
q_depart = sol_chief(t_dep_start)

# Apply change in velocity at that point
dv = @SVector [0.0, 0.0, 0.0, 0.00017947635829913735, -5.932764319639683e-5, -0.0009575628577891083]
q_deputy = q_depart + dv

# Build the callback
threshold = 0.00026014568158168577
cb = distance_callback(sol_chief, threshold)

# Build & Solve the deputy problem
tspan_deputy = (t_dep_start, tof_chief - t_dep_start)
prob_deputy = ODEProblem{false}(eoms, q_deputy, tspan_deputy, μ)
sol_deputy = solve(prob_deputy, Vern9(); callback=cb, abstol=1e-12, reltol=1e-12)

gravity = 9.8
stiffness = 500
equilibrium_length = 1
T = 5.0

f2(u, p, t) = begin
    x1, x2, dx1, dx2 = u
    length = abs(x2 - x1)
    spring_force = stiffness * (equilibrium_length - length)
    ddx1 = -gravity - spring_force
    if x1 <= 0
      ddx1 = max(0, ddx1)
    end
    ddx2 = -gravity + spring_force
    [dx1, dx2, ddx1, ddx2]
end

# https://github.com/SciML/DifferentialEquations.jl/issues/647
sol = solve(
    ODEProblem(f, [5.0, 6.0, 0.0, 0.0], (0.0, T)),
    Tsit5(),
    callback = ContinuousCallback((u, _, _) -> u[1], (integrator) -> (integrator.u[3] = 0)),
    reltol = 1e-2,
    abstol = 1e-2
)

# https://github.com/SciML/DifferentialEquations.jl/issues/601
f2(u, p, t) = return [u[2], 0]
prob = ODEProblem(f, [-10., 1.], (-10.0, 10.0), nothing)

is_wall(u, t, integrator) = return u[1]
affect!(integrator) = println("Wall!")
cb = ContinuousCallback(is_wall,
                            affect!
                            )

sol = solve(prob, Tsit5(), callback = cb)

# https://github.com/SciML/DiffEqBase.jl/issues/599

# The actual system of ODEs.  This is a simple system of ODEs here.
function eom!(du,u,params,t)
    issaturated1, issaturated2 = params
    du[1] = u[2] * !Bool(issaturated1)
    du[2] = 3.0 * !any(Bool.([issaturated1, issaturated2]))
    nothing
end

# This condition function is continuously checked throughout
# the whole time-domain simulation.
function condition!(out, u, t, integrator, limits)
    if integrator.p[1] == false            #if state u1 is within saturation limits
        out[1] =  (u[1] - limits[1][2]);   #up crossing of state 1 upper bound
        out[2] = -(u[1] - limits[1][1]);   #up crossing of state 1 lower bound
    else
        out[1] = -1
        out[2] = -1
    end
end

# The "affect!" function below seems to only be called
# when either out[1] or out[2] is zero in the condition! function above.
# That is, when either the state u1 has reached a lower or upper saturation
# limit, this "affect!" function is called.  printing "event_index" will show
# you if out[1] was zero or out[2] was zero.
function affect!(integrator, event_index)
	println("\ninside affect! function..");
    @show event_index
    if event_index ∈ [1,2] # 1 of the 2 sat. limits was met.
        # Set velocity to zero and unsaturate state 2, assumes 0 ∈ bounds of state 2
        println("integrator: ", integrator);
        println("full_cache(integrator): ", full_cache(integrator));
        println("typeof(full_cache(integrator)): ", typeof(full_cache(integrator)));
        println("integrator.p: ", integrator.p);
        println("integrator.u: ", integrator.u);
        println("integrator.du: ", integrator.du);
        println("integrator.t: ", integrator.t);
        #println("size(full_cache(integrator)): ", size(full_cache(integrator)));
        # Set velocity to zero and unsaturate state 2, assumes 0 ∈ bounds of state 2
        for u_ in full_cache(integrator)
            #println("u_: ", u_);
            u_[2] = 0.0;
            #println("u_: ", u_);
        end
        #println("u_: ", u_);
        println("full_cache(integrator): ", full_cache(integrator));
        println("integrator.p: ", integrator.p);
        #integrator.p[2] = false;
        println("integrator.p: ", integrator.p);

        u2 = deepcopy(integrator.u);
        println("u2: ", u2);
        p2 = deepcopy(integrator.p);
        p2[1:2] .= false
        du = similar(u2)
        integrator.f(du,u2,p2,integrator); # this evaluate the integrator ???

        if event_index == 2 # u1 reaches lower limit.
            integrator.p[1] = du[2] > 0.0 ? false : true
        elseif event_index == 1 # u1 reaches upper limit.
            integrator.p[1] = du[2] < 0.0 ? false : true
        end

        println("integrator.p: ", integrator.p);

    end
end

function affect2!(integrator, event_index)
    println("\ninside affect! function..");
    @show event_index
    if event_index ∈ [1,2] # 1 of the 2 sat. limits was met.
        # Set velocity to zero and unsaturate state 2, assumes 0 ∈ bounds of state 2
        println("integrator: ", integrator);
        println("full_cache(integrator): ", full_cache(integrator));
        println("typeof(full_cache(integrator)): ", typeof(full_cache(integrator)));
        println("integrator.p: ", integrator.p);
        println("integrator.u: ", integrator.u);
        println("integrator.du: ", integrator.du);
        println("integrator.t: ", integrator.t);

        integrator.p[2] = false;
        println("integrator.p: ", integrator.p);

        u2 = deepcopy(integrator.u)
        u2[2] = 0.0;
        p2 = deepcopy(integrator.p)
        p2[1:2] .= false
        du = similar(u2)
        integrator.f(du,u2,p2,integrator); # this evaluate the integrator ???

        if event_index == 2 # u1 reaches lower limit.
            integrator.p[1] = du[2] > 0.0 ? false : true
        elseif event_index == 1 # u1 reaches upper limit.
            integrator.p[1] = du[2] < 0.0 ? false : true
        end

        #integrator.p[1] = true;

        println("integrator.p: ", integrator.p);

    end
end

# ----------------------------------------
# The actual simulation is below.
# ----------------------------------------
time_to_sim = 8.0;
u0          = [0.0,-10.0];
tspan       = (0.0,time_to_sim);
data_times  = collect(0.0:1.0e-3:time_to_sim);

# Set saturation limits, compute if IC in limits, and set params accordingly
saturationlimits = [[-10.0,10.0],[-Inf,Inf]];
issaturated = .!([lim[1] <= u_ <= lim[2] for (u_,lim) ∈ zip(u0, saturationlimits)])
params = [issaturated...];

#cb = VectorContinuousCallback((out, u, t, integrator) -> condition!(out, u, t, integrator, saturationlimits),
#                              affect2!,2, affect_neg! = nothing, save_positions=(true,true));

cb = VectorContinuousCallback((out, u, t, integrator) -> condition!(out, u, t, integrator, saturationlimits),
                            affect2!,2, save_positions=(true,true));

prob_sat = ODEProblem(eom!,u0,tspan,params);
sol_nocb = solve(prob_sat,Tsit5(), saveat=data_times);
sol_cb   = solve(prob_sat,Tsit5(), callback = cb, saveat=data_times);
