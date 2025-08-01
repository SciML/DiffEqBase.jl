module DiffEqBaseTrackerExt

using DiffEqBase
import DiffEqBase: value
import Tracker

DiffEqBase.value(x::Type{Tracker.TrackedReal{T}}) where {T} = T
DiffEqBase.value(x::Type{Tracker.TrackedArray{T, N, A}}) where {T, N, A} = Array{T, N}
DiffEqBase.value(x::Tracker.TrackedReal) = x.data
DiffEqBase.value(x::Tracker.TrackedArray) = x.data

DiffEqBase.unitfulvalue(x::Type{Tracker.TrackedReal{T}}) where {T} = T
function DiffEqBase.unitfulvalue(x::Type{Tracker.TrackedArray{T, N, A}}) where {T, N, A}
    Array{T, N}
end
DiffEqBase.unitfulvalue(x::Tracker.TrackedReal) = x.data
DiffEqBase.unitfulvalue(x::Tracker.TrackedArray) = x.data

DiffEqBase.promote_u0(u0::Tracker.TrackedArray, p::Tracker.TrackedArray, t0) = u0
function DiffEqBase.promote_u0(u0::AbstractArray{<:Tracker.TrackedReal},
        p::Tracker.TrackedArray, t0)
    u0
end
function DiffEqBase.promote_u0(u0::Tracker.TrackedArray,
        p::AbstractArray{<:Tracker.TrackedReal}, t0)
    u0
end
function DiffEqBase.promote_u0(u0::AbstractArray{<:Tracker.TrackedReal},
        p::AbstractArray{<:Tracker.TrackedReal}, t0)
    u0
end
DiffEqBase.promote_u0(u0, p::Tracker.TrackedArray, t0) = Tracker.track(u0)
DiffEqBase.promote_u0(u0, p::AbstractArray{<:Tracker.TrackedReal}, t0) = eltype(p).(u0)

@inline Base.any(f::Function, x::Tracker.TrackedArray) = any(f, Tracker.data(x))

# Support adaptive with non-tracked time
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::Tracker.TrackedArray, t)
    sqrt(sum(abs2, DiffEqBase.value(u)) / length(u))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::AbstractArray{<:Tracker.TrackedReal, N},
        t) where {N}
    sqrt(sum(x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]),
        zip((DiffEqBase.value(x) for x in u), Iterators.repeated(t))) / length(u))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::Array{<:Tracker.TrackedReal, N},
        t) where {N}
    sqrt(sum(x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]),
        zip((DiffEqBase.value(x) for x in u), Iterators.repeated(t))) / length(u))
end
@inline DiffEqBase.ODE_DEFAULT_NORM(u::Tracker.TrackedReal, t) = abs(DiffEqBase.value(u))

# Support TrackedReal time, don't drop tracking on the adaptivity there
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::Tracker.TrackedArray,
        t::Tracker.TrackedReal)
    sqrt(sum(abs2, u) / length(u))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::AbstractArray{<:Tracker.TrackedReal, N},
        t::Tracker.TrackedReal) where {N}
    sqrt(sum(x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]), zip(u, Iterators.repeated(t))) /
         length(u))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::Array{<:Tracker.TrackedReal, N},
        t::Tracker.TrackedReal) where {N}
    sqrt(sum(x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]), zip(u, Iterators.repeated(t))) /
         length(u))
end
@inline DiffEqBase.ODE_DEFAULT_NORM(u::Tracker.TrackedReal, t::Tracker.TrackedReal) = abs(u)

function DiffEqBase.solve_up(prob::DiffEqBase.AbstractDEProblem,
        sensealg::Union{
            SciMLBase.AbstractOverloadingSensitivityAlgorithm,
            Nothing}, u0::Tracker.TrackedArray,
        p::Tracker.TrackedArray, args...; kwargs...)
    Tracker.track(DiffEqBase.solve_up, prob, sensealg, u0, p, args...; kwargs...)
end

function DiffEqBase.solve_up(prob::DiffEqBase.AbstractDEProblem,
        sensealg::Union{
            SciMLBase.AbstractOverloadingSensitivityAlgorithm,
            Nothing}, u0::Tracker.TrackedArray, p, args...;
        kwargs...)
    Tracker.track(DiffEqBase.solve_up, prob, sensealg, u0, p, args...; kwargs...)
end

function DiffEqBase.solve_up(prob::DiffEqBase.AbstractDEProblem,
        sensealg::Union{
            SciMLBase.AbstractOverloadingSensitivityAlgorithm,
            Nothing}, u0, p::Tracker.TrackedArray, args...;
        kwargs...)
    Tracker.track(DiffEqBase.solve_up, prob, sensealg, u0, p, args...; kwargs...)
end

Tracker.@grad function DiffEqBase.solve_up(prob,
        sensealg::Union{Nothing,
            SciMLBase.AbstractOverloadingSensitivityAlgorithm
        },
        u0, p, args...;
        kwargs...)
    sol,
    pb_f = DiffEqBase._solve_adjoint(
        prob, sensealg, Tracker.data(u0), Tracker.data(p),
        SciMLBase.TrackerOriginator(), args...; kwargs...)

    if sol isa AbstractArray
        !hasfield(typeof(sol), :u) && return sol, pb_f # being safe here
        return sol.u, pb_f # AbstractNoTimeSolution isa AbstractArray
    end
    return convert(AbstractArray, sol), pb_f
end

end
