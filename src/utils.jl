# Handled in Extensions
value(x) = x
isdistribution(u0) = false

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

# `reduce` and `map` are specialized on tuples to be unrolled (via recursion)
# Therefore, they can be type stable even with heterogeneous input types.
# We also don't care about allocating any temporaries with them, as it should
# all be unrolled and optimized away.
# Being unrolled also means const prop can work for things like
# `mapreduce(f, op, propertynames(x))`
# where `f` may call `getproperty` and thus have return type dependent
# on the particular symbol.
# `mapreduce` hasn't received any such specialization.
@inline diffeqmapreduce(f::F, op::OP, x::Tuple) where {F, OP} = reduce_tup(op, map(f, x))
@inline function diffeqmapreduce(f::F, op::OP, x::NamedTuple) where {F, OP}
    reduce_tup(op, map(f, x))
end
# For other container types, we probably just want to call `mapreduce`
@inline diffeqmapreduce(f::F, op::OP, x) where {F, OP} = mapreduce(f, op, x, init = Any)

struct DualEltypeChecker{T, T2}
    x::T
    counter::T2
end

getval(::Val{I}) where {I} = I
getval(::Type{Val{I}}) where {I} = I
getval(I::Int) = I

function (dec::DualEltypeChecker)(::Val{Y}) where {Y}
    isdefined(dec.x, Y) || return Any
    getval(dec.counter) >= DUALCHECK_RECURSION_MAX && return Any
    anyeltypedual(getfield(dec.x, Y), Val{getval(dec.counter)})
end

# Untyped dispatch: catch composite types, check all of their fields
"""
    anyeltypedual(x)


Searches through a type to see if any of its values are parameters. This is used to
then promote other values to match the dual type. For example, if a user passes a parameter

which is a `Dual` and a `u0` which is a `Float64`, after the first time step, `f(u,p,t) = p*u`
will change `u0` from `Float64` to `Dual`. Thus the state variable always needs to be converted
to a dual number before the solve. Worse still, this needs to be done in the case of
`f(du,u,p,t) = du[1] = p*u[1]`, and thus running `f` and taking the return value is not a valid
way to calculate the required state type.

But given the properties of automatic differentiation requiring that differentiation of parameters
implies differentiation of state, we assume any dual parameters implies differentiation of state
and then attempt to upconvert `u0` to match that dual-ness. Because this changes types, this needs
to be specified at compiled time and thus cannot have a Bool-based opt out, so in the future this
may be extended to use a preference system to opt-out with a `UPCONVERT_DUALS`. In the case where
upconversion is not done automatically, the user is required to upconvert all initial conditions
themselves, for an example of how this can be confusing to a user see
https://discourse.julialang.org/t/typeerror-in-julia-turing-when-sampling-for-a-forced-differential-equation/82937
"""
@generated function anyeltypedual(x, ::Type{Val{counter}} = Val{0}) where {counter}
    x = x.name === Core.Compiler.typename(Type) ? x.parameters[1] : x
    if x <: ForwardDiff.Dual
        :($x)
    elseif fieldnames(x) === ()
        :(Any)
    elseif counter < DUALCHECK_RECURSION_MAX
        T = diffeqmapreduce(x -> anyeltypedual(x, Val{counter + 1}), promote_dual,
            x.parameters)
        if T === Any || isconcretetype(T)
            :($T)
        else
            :(diffeqmapreduce(DualEltypeChecker($x, $counter + 1), promote_dual,
                map(Val, fieldnames($(typeof(x))))))
        end
    else
        :(Any)
    end
end

const FORWARDDIFF_AUTODETECTION_FAILURE_MESSAGE = """
                                Failed to automatically detect ForwardDiff compatability of
                                the parameter object. In order for ForwardDiff.jl automatic
                                differentiation to work on a solution object, the state of
                                the differential equation or nonlinear solve (`u0`) needs to
                                be converted to a Dual type which matches the values being
                                differentiated. For example, for a loss function loss(p)
                                where `p`` is a `Vector{Float64}`, this conversion is
                                equivalent to:

                                ```julia
                                # Convert u0 to match the new Dual element type of `p`
                                _prob = remake(prob, u0 = eltype(p).(prob.u0))
                                ```

                                In most cases, SciML tools are able to do this conversion
                                automatically. However, it seems you have provided a
                                parameter type for which this automatic conversion has failed.

                                To fix this, you can do the conversion yourself. For example,
                                if you have a parameter vector being optimized `p` which is
                                then put into an odd struct, you can manually convert `u0`
                                to match `p`:

                                ```julia
                                function loss(p)
                                    _prob = remake(prob, u0 = eltype(p).(prob.u0), p = MyStruct(p))
                                    sol = solve(_prob, ...)
                                    # do stuff on sol
                                end
                                ```

                                Or you can define a dispatch on `DiffEqBase.anyeltypedual`
                                which tells the system what fields to interpret as the
                                differentiable parts. For example, to support ODESolutions
                                as parameters we tell it the data is `sol.u` and `sol.t` via:

                                ```julia
                                function DiffEqBase.anyeltypedual(sol::ODESolution, counter = 0)
                                    DiffEqBase.anyeltypedual((sol.u, sol.t))
                                end
                                ```

                                To opt a type out of the dual checking, define an overload
                                that returns Any. For example:

                                ```julia
                                function DiffEqBase.anyeltypedual(::YourType, ::Type{Val{counter}}) where {counter}
                                    Any
                                end
                                ```

                                If you have defined this on a common type which should
                                be more generally supported, please open a pull request
                                adding this dispatch. If you need help defining this dispatch,
                                feel free to open an issue.
                                """

struct ForwardDiffAutomaticDetectionFailure <: Exception end

function Base.showerror(io::IO, e::ForwardDiffAutomaticDetectionFailure)
    print(io, FORWARDDIFF_AUTODETECTION_FAILURE_MESSAGE)
end

function anyeltypedual(::Type{Union{}})
    throw(ForwardDiffAutomaticDetectionFailure())
end

function anyeltypedual(::Type{<:AbstractTimeseriesSolution{T, N}},
        ::Type{Val{counter}} = Val{0}) where {T, N, counter}
    anyeltypedual(T)
end

function anyeltypedual(
        ::Type{T},
        ::Type{Val{counter}} = Val{0}) where {counter} where {T <:
                                                              NonlinearProblem{
        uType, iip, pType}} where {uType, iip, pType}
    return anyeltypedual((uType, pType), Val{counter})
end

function anyeltypedual(
        ::Type{T},
        ::Type{Val{counter}} = Val{0}) where {counter} where {T <:
                                                              NonlinearLeastSquaresProblem{
        uType, iip, pType}} where {uType, iip, pType}
    return anyeltypedual((uType, pType), Val{counter})
end

function anyeltypedual(x::SciMLBase.RecipesBase.AbstractPlot,
        ::Type{Val{counter}} = Val{0}) where {counter}
    Any
end
function anyeltypedual(x::Returns, ::Type{Val{counter}} = Val{0}) where {counter}
    anyeltypedual(x.value, Val{counter})
end

Base.@assume_effects :foldable function __anyeltypedual(::Type{T}) where {T}
    if T isa Union
        promote_dual(anyeltypedual(T.a), anyeltypedual(T.b))
    elseif hasproperty(T, :parameters)
        mapreduce(anyeltypedual, promote_dual, T.parameters; init = Any)
    else
        T
    end
end
function anyeltypedual(::Type{T}, ::Type{Val{counter}} = Val{0}) where {counter} where {T}
    __anyeltypedual(T)
end

function anyeltypedual(::Type{T},
    ::Type{Val{counter}} = Val{0}) where {counter} where {T <:
                                                          Union{AbstractArray, Set}}
anyeltypedual(eltype(T))
end
Base.@pure function __anyeltypedual_ntuple(::Type{T}) where {T <: NTuple}
if isconcretetype(eltype(T))
    return eltype(T)
end
if isempty(T.parameters)
    Any
else
    mapreduce(anyeltypedual, promote_dual, T.parameters; init = Any)
end
end
function anyeltypedual(
    ::Type{T}, ::Type{Val{counter}} = Val{0}) where {counter} where {T <: NTuple}
__anyeltypedual_ntuple(T)
end

# Any in this context just means not Dual
function anyeltypedual(
    x::SciMLBase.NullParameters, ::Type{Val{counter}} = Val{0}) where {counter}
Any
end

function anyeltypedual(sol::RecursiveArrayTools.AbstractDiffEqArray, counter = 0)
diffeqmapreduce(anyeltypedual, promote_dual, (sol.u, sol.t))
end

function anyeltypedual(prob::Union{ODEProblem, SDEProblem, RODEProblem, DDEProblem},
    ::Type{Val{counter}} = Val{0}) where {counter}
anyeltypedual((prob.u0, prob.p, prob.tspan))
end

function anyeltypedual(
    prob::Union{NonlinearProblem, NonlinearLeastSquaresProblem, OptimizationProblem},
    ::Type{Val{counter}} = Val{0}) where {counter}
anyeltypedual((prob.u0, prob.p))
end

function anyeltypedual(x::Number, ::Type{Val{counter}} = Val{0}) where {counter}
anyeltypedual(typeof(x))
end
function anyeltypedual(
    x::Union{String, Symbol}, ::Type{Val{counter}} = Val{0}) where {counter}
typeof(x)
end
function anyeltypedual(x::Union{AbstractArray{T}, Set{T}},
    ::Type{Val{counter}} = Val{0}) where {counter} where {
    T <:
    Union{Number,
    Symbol,
    String}}
anyeltypedual(T)
end
function anyeltypedual(x::Union{AbstractArray{T}, Set{T}},
    ::Type{Val{counter}} = Val{0}) where {counter} where {
    T <: Union{
    AbstractArray{
        <:Number,
    },
    Set{
        <:Number,
    }}}
anyeltypedual(eltype(x))
end
function anyeltypedual(x::Union{AbstractArray{T}, Set{T}},
    ::Type{Val{counter}} = Val{0}) where {counter} where {N, T <: NTuple{N, <:Number}}
anyeltypedual(eltype(x))
end

# Try to avoid this dispatch because it can lead to type inference issues when !isconcrete(eltype(x))
function anyeltypedual(x::AbstractArray, ::Type{Val{counter}} = Val{0}) where {counter}
if isconcretetype(eltype(x))
    anyeltypedual(eltype(x))
elseif !isempty(x) && all(i -> isassigned(x, i), 1:length(x)) &&
       counter < DUALCHECK_RECURSION_MAX
    _counter = Val{counter + 1}
    mapreduce(y -> anyeltypedual(y, _counter), promote_dual, x)
else
    # This fallback to Any is required since otherwise we cannot handle `undef` in all cases
    #  misses cases of
    Any
end
end

function anyeltypedual(x::Set, ::Type{Val{counter}} = Val{0}) where {counter}
if isconcretetype(eltype(x))
    anyeltypedual(eltype(x))
else
    # This fallback to Any is required since otherwise we cannot handle `undef` in all cases
    Any
end
end

function anyeltypedual(x::Tuple, ::Type{Val{counter}} = Val{0}) where {counter}
# Handle the empty tuple case separately for inference and to avoid mapreduce error
if x === ()
    Any
else
    diffeqmapreduce(anyeltypedual, promote_dual, x)
end
end
function anyeltypedual(x::Dict, ::Type{Val{counter}} = Val{0}) where {counter}
isempty(x) ? eltype(values(x)) : mapreduce(anyeltypedual, promote_dual, values(x))
end
function anyeltypedual(x::NamedTuple, ::Type{Val{counter}} = Val{0}) where {counter}
anyeltypedual(values(x))
end

function anyeltypedual(
    f::SciMLBase.AbstractSciMLFunction, ::Type{Val{counter}}) where {counter}
Any
end

anyeltypedual(::@Kwargs{}, ::Type{Val{counter}} = Val{0}) where {counter} = Any
anyeltypedual(::Type{@Kwargs{}}, ::Type{Val{counter}} = Val{0}) where {counter} = Any

@inline function promote_u0(u0, p, t0)
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
    Tp = anyeltypedual(p)
    if Tp == Any
        Tp = Tu
    end
    Tt = anyeltypedual(t0)
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

@inline function promote_u0(u0::AbstractArray{<:Complex}, p, t0)
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
    Tp = anyeltypedual(p)
    if Tp == Any
        Tp = Tu
    end
    Tt = anyeltypedual(t0)
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
