module DiffEqBaseForwardDiffExt

using DiffEqBase, ForwardDiff
using SimpleNonlinearSolve: ITP
using DiffEqBase.ArrayInterface
using DiffEqBase: Void, FunctionWrappersWrappers, OrdinaryDiffEqTag,
                  AbstractTimeseriesSolution,
                  RecursiveArrayTools, reduce_tup, _promote_tspan, has_continuous_callback
import DiffEqBase: hasdualpromote, wrapfun_oop, wrapfun_iip, prob2dtmin,
                   promote_tspan, ODE_DEFAULT_NORM
import SciMLBase: isdualtype, DualEltypeChecker, sse, __sum

const dualT = ForwardDiff.Dual{ForwardDiff.Tag{OrdinaryDiffEqTag, Float64}, Float64, 1}
dualgen(::Type{T}) where {T} = ForwardDiff.Dual{ForwardDiff.Tag{OrdinaryDiffEqTag, T}, T, 1}

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

function hasdualpromote(u0, t::Number)
    hasmethod(ArrayInterface.promote_eltype,
            Tuple{Type{typeof(u0)}, Type{dualgen(eltype(u0))}}) &&
        hasmethod(promote_rule,
            Tuple{Type{eltype(u0)}, Type{dualgen(eltype(u0))}}) &&
        hasmethod(promote_rule,
            Tuple{Type{eltype(u0)}, Type{typeof(t)}})
end

@inline ODE_DEFAULT_NORM(u::ForwardDiff.Dual, ::Any) = sqrt(sse(u))
@inline function ODE_DEFAULT_NORM(u::AbstractArray{<:ForwardDiff.Dual{Tag, T}},
        t::Any) where {Tag, T}
    sqrt(DiffEqBase.__sum(sse, u; init = sse(zero(T))) / DiffEqBase.totallength(u))
end
@inline ODE_DEFAULT_NORM(u::ForwardDiff.Dual, ::ForwardDiff.Dual) = sqrt(sse(u))
@inline function ODE_DEFAULT_NORM(u::AbstractArray{<:ForwardDiff.Dual{Tag, T}},
        ::ForwardDiff.Dual) where {Tag, T}
    sqrt(DiffEqBase.__sum(sse, u; init = sse(zero(T))) / DiffEqBase.totallength(u))
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

# Differentiation of internal solver

function scalar_nlsolve_ad(prob, alg::ITP, args...; kwargs...)
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
        alg::ITP, args...;
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
        alg::ITP, args...;
        kwargs...) where {uType, iip, T, V, P}
    sol, partials = scalar_nlsolve_ad(prob, alg, args...; kwargs...)

    return SciMLBase.build_solution(prob, alg, ForwardDiff.Dual{T, V, P}(sol.u, partials),
        sol.resid; retcode = sol.retcode,
        left = ForwardDiff.Dual{T, V, P}(sol.left, partials),
        right = ForwardDiff.Dual{T, V, P}(sol.right, partials))
end

end
