module DiffEqBaseMonteCarloMeasurementsExt

if isdefined(Base, :get_extension)
    using DiffEqBase
    import DiffEqBase: value
    using MonteCarloMeasurements
else
    using ..DiffEqBase
    import ..DiffEqBase: value
    using ..MonteCarloMeasurements
end

function DiffEqBase.promote_u0(u0::AbstractArray{
        <:MonteCarloMeasurements.AbstractParticles,
    },
    p::AbstractArray{<:MonteCarloMeasurements.AbstractParticles},
    t0)
    u0
end
function DiffEqBase.promote_u0(u0,
    p::AbstractArray{<:MonteCarloMeasurements.AbstractParticles},
    t0)
    eltype(p).(u0)
end

DiffEqBase.value(x::Type{MonteCarloMeasurements.AbstractParticles{T, N}}) where {T, N} = T
DiffEqBase.value(x::MonteCarloMeasurements.AbstractParticles) = mean(x.particles)

@inline function DiffEqBase.fastpow(x::MonteCarloMeasurements.AbstractParticles,
    y::MonteCarloMeasurements.AbstractParticles)
    x^y
end

# Support adaptive steps should be errorless
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::AbstractArray{
        <:MonteCarloMeasurements.AbstractParticles,
        N}, t) where {N}
    sqrt(mean(x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]),
        zip((value(x) for x in u), Iterators.repeated(t))))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::AbstractArray{
        <:MonteCarloMeasurements.AbstractParticles,
        N},
    t::AbstractArray{
        <:MonteCarloMeasurements.AbstractParticles,
        N}) where {N}
    sqrt(mean(x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]),
        zip((value(x) for x in u), Iterators.repeated(value.(t)))))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::MonteCarloMeasurements.AbstractParticles, t)
    abs(value(u))
end

# To make DAE problems work, we translate the MCM problem into an ensemble problem,
# solve it, and then repackage the solution into an ODESolution with Particles

const ParticleODEProblem = ODEProblem{<:Any, <:Any, <:Any, <:MonteCarloMeasurements.SomeKindOfParticles}

using MonteCarloMeasurements: nparticles, vecindex

function SciMLBase.EnsembleProblem(prob::ParticleODEProblem, args...; kwargs...)
    p = copy(prob.p)
    u0 = copy(prob.u0)
    N = max(nparticles(eltype(p)), nparticles(eltype(u0)))
    p0 = vecindex.(p, 1)
    u00 = vecindex.(u0, 1)
    prob0 = remake(prob, p = p0, u0 = u00)

    prob_func = function(probi,i,repeat)
        for j in eachindex(probi.u0)
            probi.u0[j] = vecindex(u0[j], i)
        end
        for j in eachindex(probi.p)
            probi.p[j] = vecindex(p[j], i)
        end
        probi
    end
    
    eprob = EnsembleProblem(prob0; prob_func)
end

function DiffEqBase.solve(prob::ParticleODEProblem, alg, args...; kwargs...)
    N = max(nparticles(eltype(prob.p)), nparticles(eltype(prob.u0)))
    eprob = EnsembleProblem(prob)
    esol = DiffEqBase.solve(eprob, alg, EnsembleThreads(); trajectories = N, kwargs...)
    postprocess(esol)
end

function postprocess(esol)
    soli = esol[1]
    t = soli.t
    nx = length(soli[1])

    # [state_index][mc_index]
    utraj = map(t) do t
        data = reduce(hcat, OrdinaryDiffEq.EnsembleAnalysis.componentwise_vectors_timepoint(esol,t)) # nmc × nx
        xparts = Particles(data)
    end

    u_analytic  = nothing
    errors      = nothing
    k           = nothing
    prob        = soli.prob
    alg         = soli.alg
    interp      = nothing
    dense       = false
    tslocation  = soli.tslocation
    destats     = nothing
    alg_choice  = soli.alg_choice
    retcode     = soli.retcode
    ODESolution{eltype(soli.u[1]), 2}(utraj, u_analytic, errors, t, k, prob, alg, interp, dense,
    tslocation, destats, alg_choice, retcode)
end

end
