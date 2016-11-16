# This is a simple example showing how forward and backward Euler
# could be wrapped

using DiffEqBase
const D=DiffEqBase

# make a algorithm type
abstract EulerAlgs <: AbstractODEAlgorithm

immutable FwdEulerAlg <: EulerAlgs
end
# register it
D.add_alg!(FwdEulerAlg, "Forward Euler", D.notgreat, D.nonstiff, D.low)

immutable BwdEulerAlg <: EulerAlgs
end
# register it
D.add_alg!(BwdEulerAlg, "Backward Euler", D.notgreat, D.stiff, D.low)

function D.solve{uType,tType,isinplace}(p::AbstractODEProblem{uType,tType,isinplace},
                                                 Alg::Type{FwdEulerAlg};
                                                 dt=(p.tspan[2]-p.tspan[1])/100,
                                                 tstops=tType[],
                                                 kwargs... # ignored kwargs
                                                 )
    u0 = p.u0
    f = p.f
    tspan = p.tspan
    @assert !isinplace "Only out of place functions supported"

    if isempty(tstops)
        tstops = tspan[1]:dt:tspan[2]
    end
    @assert tstops[1]==tspan[1]

    nt = length(tstops)
    out = Vector{uType}(nt)
    out[1] = copy(u0)
    for i=2:nt
        t = tstops[i]
        dt = t-tstops[i-1]
        out[i] = out[i-1] + dt*f(t,out[i-1])
    end
    # make solution type
    build_ode_solution(p, Alg(), tstops, out)
end


# Try it
dt = 0.01
p1 = ODETestProblem((t,u)->u, 1.0, (t,u) -> exp(t), (0.0,1.0))
sol1 = solve(p1, FwdEulerAlg, tstops=0:dt:1)

function D.solve{uType,tType,isinplace}(p::AbstractODEProblem{uType,tType,isinplace},
                                                 Alg::Type{BwdEulerAlg};
                                                 dt=(p.tspan[2]-p.tspan[1])/100,
                                                 tstops=tType[],
                                                 tol=1e-5,
                                                 maxiter=100,
                                                 kwargs... # ignored kwargs
                                                 )
    u0 = p.u0
    f = p.f
    tspan = p.tspan
    # TODO: fix numparameters as it picks up the Jacobian
#    @assert !isinplace "Only out of place functions supported"
    @assert D.has_jac(f) "Provide Jacobian as f(::Val{:jac}, ...)"
    jac = (t,u) -> f(Val{:jac}(), t, u)

    if isempty(tstops)
        tstops = tspan[1]:dt:tspan[2]
    end
    @assert tstops[1]==tspan[1]

    nt = length(tstops)
    out = Vector{uType}(nt)
    out[1] = copy(u0)
    for i=2:nt
        t = tstops[i]
        dt = t-tstops[i-1]
        out[i] = newton(t, dt, out[i-1], f, jac, tol, maxiter)
    end
    # make solution type
    build_ode_solution(p, Alg(), tstops, out)
end

function newton(t, dt, u_last, f, jac, tol, maxiter)
    res = (u) -> u - u_last - dt*f(t,u)
    u = u_last + dt*f(t,u_last) # forward Euler step as first guess
    for i=1:maxiter
        du = -(I - dt*jac(t,u))\res(u)
        u += du
        norm(du, Inf)<tol && return u
    end
    error("Newton not converged")
end

# Try it
ff(t,u) = u
ff(::Val{:jac}, t, u) = 1
p2 = ODETestProblem(ff, 1.0, (t,u) -> exp(t), (0.0,1.0) )
sol2 = solve(p2, BwdEulerAlg, tstops=0:dt:1)

### try automatic solve selection:
solve(p1, class=D.nonstiff, precision=D.low)
solve(p2, class=D.stiff, precision=D.low)

using Plots
plot(sol1)
plot!(sol2)
#plot!(sol2, plot_analytic=true)

# using DiffEqDevTools
# dts = 1./2.^(8:-1:4)
# sim = test_convergence(dts,p2,BwdEulerAlg)
# @show sim.𝒪est[:final]
# plot(sim)
