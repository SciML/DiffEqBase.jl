# This is a simple example showing how forward and backward Euler
# could be wrapped

module InternalEuler

using DiffEqBase

# make a algorithm type
abstract type EulerAlgs <: DiffEqBase.AbstractODEAlgorithm end

struct FwdEulerAlg <: EulerAlgs end

struct BwdEulerAlg <: EulerAlgs end

function DiffEqBase.solve{uType,tType,isinplace}(p::AbstractODEProblem{uType,tType,isinplace},
                                                 Alg::FwdEulerAlg;
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
    build_solution(p, Alg, tstops, out)
end

function DiffEqBase.solve{uType,tType,isinplace}(p::AbstractODEProblem{uType,tType,isinplace},
                                                 Alg::BwdEulerAlg;
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
    @assert DiffEqBase.has_jac(f) "Provide Jacobian as f(::Val{:jac}, ...)"
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
    build_solution(p, Alg, tstops, out)
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

end
