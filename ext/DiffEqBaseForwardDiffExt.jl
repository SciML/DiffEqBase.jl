module DiffEqBaseForwardDiffExt

using DiffEqBase, ForwardDiff
using DiffEqBase: Void, FunctionWrappersWrappers, OrdinaryDiffEqTag
import DiffEqBase: hasdualpromote, wrapfun_oop, wrapfun_iip, promote_u0, prob2dtmin,
promote_tspan, anyeltypedual, isdualtype, value, ODE_DEFAULT_NORM, InternalITP,
nextfloat_tdir

const DUALCHECK_RECURSION_MAX = 10

eltypedual(x)  = eltype(x) <: ForwardDiff.Dual
isdualtype(::Type{<:ForwardDiff.Dual}) = true 
const dualT = ForwardDiff.Dual{ForwardDiff.Tag{OrdinaryDiffEqTag, Float64}, Float64, 1}
dualgen(::Type{T}) where {T} = ForwardDiff.Dual{ForwardDiff.Tag{OrdinaryDiffEqTag, T}, T, 1}

# Copy of the other prob2dtmin dispatch, just for optionality
function prob2dtmin(tspan, ::ForwardDiff.Dual, use_end_time)
    t1, t2 = tspan
    isfinite(t1) || throw(ArgumentError("t0 in the tspan `(t0, t1)` must be finite"))
    if use_end_time && isfinite(t2 - t1)
        return max(eps(t2), eps(t1))
    else
        return max(eps(typeof(t1)), eps(t1))
    end
end

hasdualpromote(u0,t::Number) = hasmethod(ArrayInterface.promote_eltype,
                    Tuple{Type{typeof(u0)}, Type{dualgen(eltype(u0))}}) &&
                hasmethod(promote_rule,
                    Tuple{Type{eltype(u0)}, Type{dualgen(eltype(u0))}}) &&
                hasmethod(promote_rule,
                    Tuple{Type{eltype(u0)}, Type{typeof(t)}})

const NORECOMPILE_IIP_SUPPORTED_ARGS = (
    Tuple{Vector{Float64}, Vector{Float64},
        Vector{Float64}, Float64},
    Tuple{Vector{Float64}, Vector{Float64},
        SciMLBase.NullParameters, Float64})

const oop_arglists = (Tuple{Vector{Float64}, Vector{Float64}, Float64},
    Tuple{Vector{Float64}, SciMLBase.NullParameters, Float64},
    Tuple{Vector{Float64}, Vector{Float64}, dualT},
    Tuple{Vector{dualT}, Vector{Float64}, Float64},
    Tuple{Vector{dualT}, SciMLBase.NullParameters, Float64},
    Tuple{Vector{Float64}, SciMLBase.NullParameters, dualT})

const NORECOMPILE_OOP_SUPPORTED_ARGS = (Tuple{Vector{Float64},
        Vector{Float64}, Float64},
    Tuple{Vector{Float64},
        SciMLBase.NullParameters, Float64})
const oop_returnlists = (Vector{Float64}, Vector{Float64},
    ntuple(x -> Vector{dualT}, length(oop_arglists) - 2)...)

function wrapfun_oop(ff, inputs::Tuple = ())
    if !isempty(inputs)
        IT = Tuple{map(typeof, inputs)...}
        if IT ∉ NORECOMPILE_OOP_SUPPORTED_ARGS
            throw(NoRecompileArgumentError(IT))
        end
    end
    FunctionWrappersWrappers.FunctionWrappersWrapper(ff, oop_arglists,
        oop_returnlists)
end

function wrapfun_iip(ff,
        inputs::Tuple{T1, T2, T3, T4}) where {T1, T2, T3, T4}
    T = eltype(T2)
    dualT = dualgen(T)
    dualT1 = ArrayInterface.promote_eltype(T1, dualT)
    dualT2 = ArrayInterface.promote_eltype(T2, dualT)
    dualT4 = dualgen(promote_type(T, T4))

    iip_arglists = (Tuple{T1, T2, T3, T4},
        Tuple{dualT1, dualT2, T3, T4},
        Tuple{dualT1, T2, T3, dualT4},
        Tuple{dualT1, dualT2, T3, dualT4})

    iip_returnlists = ntuple(x -> Nothing, 4)

    fwt = map(iip_arglists, iip_returnlists) do A, R
        FunctionWrappersWrappers.FunctionWrappers.FunctionWrapper{R, A}(Void(ff))
    end
    FunctionWrappersWrappers.FunctionWrappersWrapper{typeof(fwt), false}(fwt)
end

const iip_arglists_default = (
    Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64},
        Float64},
    Tuple{Vector{Float64}, Vector{Float64},
        SciMLBase.NullParameters,
        Float64
    },
    Tuple{Vector{dualT}, Vector{Float64}, Vector{Float64}, dualT},
    Tuple{Vector{dualT}, Vector{dualT}, Vector{Float64}, dualT},
    Tuple{Vector{dualT}, Vector{dualT}, Vector{Float64}, Float64},
    Tuple{Vector{dualT}, Vector{dualT}, SciMLBase.NullParameters,
        Float64
    },
    Tuple{Vector{dualT}, Vector{Float64},
        SciMLBase.NullParameters, dualT
    })
const iip_returnlists_default = ntuple(x -> Nothing, length(iip_arglists_default))

function wrapfun_iip(@nospecialize(ff))
    fwt = map(iip_arglists_default, iip_returnlists_default) do A, R
        FunctionWrappersWrappers.FunctionWrappers.FunctionWrapper{R, A}(Void(ff))
    end
    FunctionWrappersWrappers.FunctionWrappersWrapper{typeof(fwt), false}(fwt)
end

"""
    promote_dual(::Type{T},::Type{T2})


Is like the number promotion system, but always prefers a dual number type above
anything else. For higher order differentiation, it returns the most dualiest of
them all. This is then used to promote `u0` into the suspected highest differentiation
space for solving the equation.
"""
promote_dual(::Type{T}, ::Type{T2}) where {T, T2} = T
promote_dual(::Type{T}, ::Type{T2}) where {T <: ForwardDiff.Dual, T2} = T
function promote_dual(::Type{T},
        ::Type{T2}) where {T <: ForwardDiff.Dual, T2 <: ForwardDiff.Dual}
    T
end
promote_dual(::Type{T}, ::Type{T2}) where {T, T2 <: ForwardDiff.Dual} = T2

function promote_dual(::Type{T},
        ::Type{T2}) where {T3, T4, V, V2 <: ForwardDiff.Dual, N, N2,
        T <: ForwardDiff.Dual{T3, V, N},
        T2 <: ForwardDiff.Dual{T4, V2, N2}}
    T2
end
function promote_dual(::Type{T},
        ::Type{T2}) where {T3, T4, V <: ForwardDiff.Dual, V2, N, N2,
        T <: ForwardDiff.Dual{T3, V, N},
        T2 <: ForwardDiff.Dual{T4, V2, N2}}
    T
end
function promote_dual(::Type{T},
        ::Type{T2}) where {
        T3, V <: ForwardDiff.Dual, V2 <: ForwardDiff.Dual,
        N,
        T <: ForwardDiff.Dual{T3, V, N},
        T2 <: ForwardDiff.Dual{T3, V2, N}}
    ForwardDiff.Dual{T3, promote_dual(V, V2), N}
end

# Opt out since these are using for preallocation, not differentiation
function anyeltypedual(x::Union{ForwardDiff.AbstractConfig, Module},
        ::Type{Val{counter}} = Val{0}) where {counter}
    Any
end
function anyeltypedual(x::Type{T},
        ::Type{Val{counter}} = Val{0}) where {counter} where {T <:
                                                              ForwardDiff.AbstractConfig}
    Any
end

function anyeltypedual(x::ForwardDiff.DiffResults.DiffResult,
        ::Type{Val{counter}} = Val{0}) where {counter}
    Any
end
function anyeltypedual(x::Type{T},
        ::Type{Val{counter}} = Val{0}) where {counter} where {T <:
                                                              ForwardDiff.DiffResults.DiffResult}
    Any
end

function anyeltypedual(::Type{T},
        ::Type{Val{counter}} = Val{0}) where {counter} where {T <: ForwardDiff.Dual}
    T
end

function promote_tspan(u0::AbstractArray{<:ForwardDiff.Dual}, p, tspan, prob, kwargs)
    if (haskey(kwargs, :callback) && has_continuous_callback(kwargs[:callback])) ||
       (haskey(prob.kwargs, :callback) && has_continuous_callback(prob.kwargs[:callback]))
        return _promote_tspan(eltype(u0).(tspan), kwargs)
    else
        return _promote_tspan(tspan, kwargs)
    end
end

function promote_tspan(u0::AbstractArray{<:Complex{<:ForwardDiff.Dual}}, p, tspan, prob,
        kwargs)
    return _promote_tspan(real(eltype(u0)).(tspan), kwargs)
end

function promote_tspan(u0::AbstractArray{<:ForwardDiff.Dual}, p,
        tspan::Tuple{<:ForwardDiff.Dual, <:ForwardDiff.Dual}, prob, kwargs)
    return _promote_tspan(tspan, kwargs)
end

value(x::Type{ForwardDiff.Dual{T, V, N}}) where {T, V, N} = V
value(x::ForwardDiff.Dual) = value(ForwardDiff.value(x))

sse(x::Number) = abs2(x)
sse(x::ForwardDiff.Dual) = sse(ForwardDiff.value(x)) + sum(sse, ForwardDiff.partials(x))
totallength(x::Number) = 1
function totallength(x::ForwardDiff.Dual)
    totallength(ForwardDiff.value(x)) + sum(totallength, ForwardDiff.partials(x))
end
totallength(x::AbstractArray) = __sum(totallength, x; init = 0)

@inline ODE_DEFAULT_NORM(u::ForwardDiff.Dual, ::Any) = sqrt(sse(u))
@inline function ODE_DEFAULT_NORM(u::AbstractArray{<:ForwardDiff.Dual{Tag, T}},
        t::Any) where {Tag, T}
    sqrt(__sum(sse, u; init = sse(zero(T))) / totallength(u))
end
@inline ODE_DEFAULT_NORM(u::ForwardDiff.Dual, ::ForwardDiff.Dual) = sqrt(sse(u))
@inline function ODE_DEFAULT_NORM(u::AbstractArray{<:ForwardDiff.Dual{Tag, T}},
        ::ForwardDiff.Dual) where {Tag, T}
    sqrt(__sum(sse, u; init = sse(zero(T))) / totallength(u))
end

if !hasmethod(nextfloat, Tuple{ForwardDiff.Dual})
    # Type piracy. Should upstream
    function Base.nextfloat(d::ForwardDiff.Dual{T, V, N}) where {T, V, N}
        ForwardDiff.Dual{T}(nextfloat(d.value), d.partials)
    end
    function Base.prevfloat(d::ForwardDiff.Dual{T, V, N}) where {T, V, N}
        ForwardDiff.Dual{T}(prevfloat(d.value), d.partials)
    end
end

# bisection(f, tup::Tuple{T,T}, t_forward::Bool) where {T<:ForwardDiff.Dual} = find_zero(f, tup, Roots.AlefeldPotraShi())

# Static Arrays don't support the `init` keyword argument for `sum`
@inline __sum(f::F, args...; init, kwargs...) where {F} = sum(f, args...; init, kwargs...)
@inline function __sum(f::F, a::DiffEqBase.StaticArraysCore.StaticArray...; init, kwargs...) where {F}
    return mapreduce(f, +, a...; init, kwargs...)
end

# Differentiation of internal solver

function scalar_nlsolve_ad(prob, alg::InternalITP, args...; kwargs...)
    f = prob.f
    p = value(prob.p)

    if prob isa IntervalNonlinearProblem
        tspan = value(prob.tspan)
        newprob = IntervalNonlinearProblem(f, tspan, p; prob.kwargs...)
    else
        u0 = value(prob.u0)
        newprob = NonlinearProblem(f, u0, p; prob.kwargs...)
    end

    sol = solve(newprob, alg, args...; kwargs...)

    uu = sol.u
    if p isa Number
        f_p = ForwardDiff.derivative(Base.Fix1(f, uu), p)
    else
        f_p = ForwardDiff.gradient(Base.Fix1(f, uu), p)
    end

    f_x = ForwardDiff.derivative(Base.Fix2(f, p), uu)
    pp = prob.p
    sumfun = let f_x′ = -f_x
        ((fp, p),) -> (fp / f_x′) * ForwardDiff.partials(p)
    end
    partials = sum(sumfun, zip(f_p, pp))
    return sol, partials
end

function SciMLBase.solve(
        prob::IntervalNonlinearProblem{uType, iip,
            <:ForwardDiff.Dual{T, V, P}},
        alg::InternalITP, args...;
        kwargs...) where {uType, iip, T, V, P}
    sol, partials = scalar_nlsolve_ad(prob, alg, args...; kwargs...)
    return SciMLBase.build_solution(prob, alg, ForwardDiff.Dual{T, V, P}(sol.u, partials),
        sol.resid; retcode = sol.retcode,
        left = ForwardDiff.Dual{T, V, P}(sol.left, partials),
        right = ForwardDiff.Dual{T, V, P}(sol.right, partials))
end

function SciMLBase.solve(
        prob::IntervalNonlinearProblem{uType, iip,
            <:AbstractArray{
                <:ForwardDiff.Dual{T,
                V,
                P},
            }},
        alg::InternalITP, args...;
        kwargs...) where {uType, iip, T, V, P}
    sol, partials = scalar_nlsolve_ad(prob, alg, args...; kwargs...)

    return SciMLBase.build_solution(prob, alg, ForwardDiff.Dual{T, V, P}(sol.u, partials),
        sol.resid; retcode = sol.retcode,
        left = ForwardDiff.Dual{T, V, P}(sol.left, partials),
        right = ForwardDiff.Dual{T, V, P}(sol.right, partials))
end

end