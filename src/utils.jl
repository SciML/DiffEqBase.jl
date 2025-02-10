# Handled in Extensions
value(x) = x
unitfulvalue(x) = x
isdistribution(u0) = false
sse(x::Number) = abs2(x)

# Static Arrays don't support the `init` keyword argument for `sum`
@inline __sum(f::F, args...; init, kwargs...) where {F} = sum(f, args...; init, kwargs...)
@inline function __sum(
        f::F, a::StaticArraysCore.StaticArray...; init, kwargs...) where {F}
    return mapreduce(f, +, a...; init, kwargs...)
end

totallength(x::Number) = 1
totallength(x::AbstractArray) = __sum(totallength, x; init = 0)

_vec(v) = vec(v)
_vec(v::Number) = v
_vec(v::AbstractSciMLScalarOperator) = v
_vec(v::AbstractVector) = v

_reshape(v, siz) = reshape(v, siz)
_reshape(v::Number, siz) = v
_reshape(v::AbstractSciMLScalarOperator, siz) = v

macro tight_loop_macros(ex)
    :($(esc(ex)))
end

# TODO: would be good to have dtmin a function of dt
function prob2dtmin(prob; use_end_time = true)
    prob2dtmin(prob.tspan, oneunit(eltype(prob.tspan)), use_end_time)
end

# This functino requires `eps` to exist, which restricts below `<: Real`
# Example of a failure is Rational
function prob2dtmin(tspan, ::AbstractFloat, use_end_time)
    t1, t2 = tspan
    isfinite(t1) || throw(ArgumentError("t0 in the tspan `(t0, t1)` must be finite"))
    if use_end_time && isfinite(t2 - t1)
        return max(eps(t2), eps(t1))
    else
        return max(eps(typeof(t1)), eps(t1))
    end
end
prob2dtmin(tspan, ::Integer, ::Any) = 0
# Multiplication is for putting the right units on the constant!
prob2dtmin(tspan, onet, ::Any) = onet * 1 // Int64(2)^33 # roughly 10^10 but more likely to turn into a multiplication.

function timedepentdtmin(integrator::DEIntegrator)
    timedepentdtmin(integrator.t, integrator.opts.dtmin)
end
timedepentdtmin(t::AbstractFloat, dtmin) = abs(max(eps(t), dtmin))
timedepentdtmin(::Any, dtmin) = abs(dtmin)

maybe_with_logger(f, logger) = logger === nothing ? f() : Logging.with_logger(f, logger)

function default_logger(logger)
    Logging.min_enabled_level(logger) ≤ ProgressLogging.ProgressLevel && return nothing

    if Sys.iswindows() || (isdefined(Main, :IJulia) && Main.IJulia.inited)
        progresslogger = ConsoleProgressMonitor.ProgressLogger()
    else
        progresslogger = TerminalLoggers.TerminalLogger()
    end

    logger1 = LoggingExtras.EarlyFilteredLogger(progresslogger) do log
        log.level == ProgressLogging.ProgressLevel
    end
    logger2 = LoggingExtras.EarlyFilteredLogger(logger) do log
        log.level != ProgressLogging.ProgressLevel
    end

    LoggingExtras.TeeLogger(logger1, logger2)
end

# for the non-unitful case the correct type is just u
_rate_prototype(u, t::T, onet::T) where {T} = u

# Nonlinear Solve functionality
@inline __fast_scalar_indexing(args...) = all(ArrayInterface.fast_scalar_indexing, args)

@inline __maximum_abs(op::F, x, y) where {F} = __maximum(abs ∘ op, x, y)
## Nonallocating version of maximum(op.(x, y))
@inline function __maximum(op::F, x, y) where {F}
    if __fast_scalar_indexing(x, y)
        return maximum(@closure((xᵢyᵢ)->begin
                xᵢ, yᵢ = xᵢyᵢ
                return op(xᵢ, yᵢ)
            end), zip(x, y))
    else
        return mapreduce(@closure((xᵢ, yᵢ)->op(xᵢ, yᵢ)), max, x, y)
    end
end

@inline function __norm_op(::typeof(Base.Fix2(norm, 2)), op::F, x, y) where {F}
    if __fast_scalar_indexing(x, y)
        return sqrt(sum(@closure((xᵢyᵢ)->begin
                xᵢ, yᵢ = xᵢyᵢ
                return op(xᵢ, yᵢ)^2
            end), zip(x, y)))
    else
        return sqrt(mapreduce(@closure((xᵢ, yᵢ)->(op(xᵢ, yᵢ)^2)), +, x, y))
    end
end

@inline __norm_op(norm::N, op::F, x, y) where {N, F} = norm(op.(x, y))

function __nonlinearsolve_is_approx(x::Number, y::Number; atol = false,
        rtol = atol > 0 ? false : sqrt(eps(promote_type(typeof(x), typeof(y)))))
    return isapprox(x, y; atol, rtol)
end
function __nonlinearsolve_is_approx(x, y; atol = false,
        rtol = atol > 0 ? false : sqrt(eps(promote_type(eltype(x), eltype(y)))))
    length(x) != length(y) && return false
    d = __maximum_abs(-, x, y)
    return d ≤ max(atol, rtol * max(maximum(abs, x), maximum(abs, y)))
end

@inline function __add_and_norm(::Nothing, x, y)
    Base.depwarn("Not specifying the internal norm of termination conditions has been \
                  deprecated. Using inf-norm currently.",
        :__add_and_norm)
    return __maximum_abs(+, x, y)
end
@inline __add_and_norm(::typeof(Base.Fix1(maximum, abs)), x, y) = __maximum_abs(+, x, y)
@inline __add_and_norm(::typeof(Base.Fix2(norm, Inf)), x, y) = __maximum_abs(+, x, y)
@inline __add_and_norm(f::F, x, y) where {F} = __norm_op(f, +, x, y)

@inline function __apply_termination_internalnorm(::Nothing, u)
    Base.depwarn("Not specifying the internal norm of termination conditions has been \
                  deprecated. Using inf-norm currently.",
        :__apply_termination_internalnorm)
    return __apply_termination_internalnorm(Base.Fix1(maximum, abs), u)
end
@inline __apply_termination_internalnorm(f::F, u) where {F} = f(u)

struct DualEltypeChecker{T, T2}
    x::T
    counter::T2
end

anyeltypedual(x) = anyeltypedual(x, Val{0})
anyeltypedual(x, counter) = Any

function promote_u0(u0, p, t0)
    if SciMLStructures.isscimlstructure(p)
        _p = SciMLStructures.canonicalize(SciMLStructures.Tunable(), p)[1]
        if !isequal(_p, p)
            return promote_u0(u0, _p, t0)
        end
    end
    Tu = eltype(u0)
    if isdualtype(Tu)
        return u0
    end
    Tp = anyeltypedual(p, Val{0})
    if Tp == Any
        Tp = Tu
    end
    Tt = anyeltypedual(t0, Val{0})
    if Tt == Any
        Tt = Tu
    end
    Tcommon = promote_type(Tu, Tp, Tt)
    return if isdualtype(Tcommon)
        Tcommon.(u0)
    else
        u0
    end
end

function promote_u0(u0::AbstractArray{<:Complex}, p, t0)
    if SciMLStructures.isscimlstructure(p)
        _p = SciMLStructures.canonicalize(SciMLStructures.Tunable(), p)[1]
        if !isequal(_p, p)
            return promote_u0(u0, _p, t0)
        end
    end
    Tu = real(eltype(u0))
    if isdualtype(Tu)
        return u0
    end
    Tp = anyeltypedual(p, Val{0})
    if Tp == Any
        Tp = Tu
    end
    Tt = anyeltypedual(t0, Val{0})
    if Tt == Any
        Tt = Tu
    end
    Tcommon = promote_type(eltype(u0), Tp, Tt)
    return if isdualtype(real(Tcommon))
        Tcommon.(u0)
    else
        u0
    end
end
