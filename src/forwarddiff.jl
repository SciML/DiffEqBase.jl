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

# `reduce` and `map` are specialized on tuples to be unrolled (via recursion)
# Therefore, they can be type stable even with heterogenous input types.
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
@inline diffeqmapreduce(f::F, op::OP, x) where {F, OP} = mapreduce(f, op, x)

struct DualEltypeChecker{T}
    x::T
    counter::Int
    DualEltypeChecker(x::T, counter::Int) where {T} = new{T}(x, counter + 1)
end
function (dec::DualEltypeChecker)(::Val{Y}) where {Y}
    isdefined(dec.x, Y) || return Any
    anyeltypedual(getproperty(dec.x, Y), dec.counter)
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

But given the properties of automatic differentiation requiring that differntiation of parameters
implies differentiation of state, we assume any dual parameters implies differentiation of state
and then attempt to upconvert `u0` to match that dual-ness. Because this changes types, this needs
to be specified at compiled time and thus cannot have a Bool-based opt out, so in the future this
may be extended to use a preference system to opt-out with a `UPCONVERT_DUALS`. In the case where
upconversion is not done automatically, the user is required to upconvert all initial conditions
themselves, for an example of how this can be confusing to a user see
https://discourse.julialang.org/t/typeerror-in-julia-turing-when-sampling-for-a-forced-differential-equation/82937
"""
function anyeltypedual(x, counter = 0)
    if propertynames(x) === ()
        Any
    elseif counter < 100
        diffeqmapreduce(DualEltypeChecker(x, counter), promote_dual,
                        map(Val, propertynames(x)))
    else
        Any
    end
end

# Opt out since these are using for preallocation, not differentiation
anyeltypedual(x::Union{ForwardDiff.AbstractConfig, Module}, counter = 0) = Any
anyeltypedual(x::Type{T}, counter = 0) where {T <: ForwardDiff.AbstractConfig} = Any
anyeltypedual(x::SciMLBase.RecipesBase.AbstractPlot, counter = 0) = Any

Base.@pure function __anyeltypedual(::Type{T}) where {T}
    hasproperty(T, :parameters) ?
    mapreduce(anyeltypedual, promote_dual, T.parameters; init = Any) : T
end
anyeltypedual(::Type{T}, counter = 0) where {T} = __anyeltypedual(T)
anyeltypedual(::Type{T}, counter = 0) where {T <: ForwardDiff.Dual} = T
function anyeltypedual(::Type{T}, counter = 0) where {T <: Union{AbstractArray, Set}}
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
anyeltypedual(::Type{T}, counter = 0) where {T <: NTuple} = __anyeltypedual_ntuple(T)

# Any in this context just means not Dual
anyeltypedual(x::SciMLBase.NullParameters, counter = 0) = Any
anyeltypedual(x::Number, counter = 0) = anyeltypedual(typeof(x))
anyeltypedual(x::Union{String, Symbol}, counter = 0) = typeof(x)
function anyeltypedual(x::Union{Array{T}, AbstractArray{T}, Set{T}},
                       counter = 0) where {
                                           T <:
                                           Union{Number,
                                                 Symbol,
                                                 String}}
    anyeltypedual(T)
end
function anyeltypedual(x::Union{Array{T}, AbstractArray{T}, Set{T}},
                       counter = 0) where {N,
                                           T <: Union{
                                                 AbstractArray{
                                                               <:Number
                                                               },
                                                 Set{
                                                     <:Number
                                                     },
                                                 NTuple{N,
                                                        <:Number
                                                        }}}
    anyeltypedual(eltype(x))
end

# Try to avoid this dispatch because it can lead to type inference issues when !isconcrete(eltype(x))
function anyeltypedual(x::AbstractArray, counter = 0)
    if isconcretetype(eltype(x))
        anyeltypedual(eltype(x))
    elseif !isempty(x) && all(i -> isassigned(x, i), 1:length(x)) && counter < 100
        counter += 1
        mapreduce(y -> anyeltypedual(y, counter), promote_dual, x)
    else
        # This fallback to Any is required since otherwise we cannot handle `undef` in all cases
        #  misses cases of
        Any
    end
end

function anyeltypedual(x::Set, counter = 0)
    if isconcretetype(eltype(x))
        anyeltypedual(eltype(x))
    else
        # This fallback to Any is required since otherwise we cannot handle `undef` in all cases
        Any
    end
end

function anyeltypedual(x::Tuple, counter = 0)
    # Handle the empty tuple case separately for inference and to avoid mapreduce error
    if x === ()
        Any
    else
        diffeqmapreduce(anyeltypedual, promote_dual, x)
    end
end
function anyeltypedual(x::Dict, counter = 0)
    isempty(x) ? eltype(values(x)) : mapreduce(anyeltypedual, promote_dual, values(x))
end
function anyeltypedual(x::NamedTuple, counter = 0)
    isempty(x) ? Any : diffeqmapreduce(anyeltypedual, promote_dual, values(x))
end
@inline function promote_u0(u0, p, t0)
    if !(eltype(u0) <: ForwardDiff.Dual)
        T = anyeltypedual(p)
        T === Any && return u0
        if T <: ForwardDiff.Dual
            return T.(u0)
        end
    end
    u0
end

function promote_tspan(u0::AbstractArray{<:ForwardDiff.Dual}, p,
                       tspan::Tuple{<:ForwardDiff.Dual, <:ForwardDiff.Dual}, prob, kwargs)
    return tspan
end

function promote_tspan(u0::AbstractArray{<:ForwardDiff.Dual}, p, tspan, prob, kwargs)
    if (haskey(kwargs, :callback) && has_continuous_callback(kwargs[:callback])) ||
       (haskey(prob.kwargs, :callback) && has_continuous_callback(prob.kwargs[:callback]))
        return eltype(u0).(tspan)
    else
        return tspan
    end
end

value(x::Type{ForwardDiff.Dual{T, V, N}}) where {T, V, N} = V
value(x::ForwardDiff.Dual) = value(ForwardDiff.value(x))

@inline fastpow(x::ForwardDiff.Dual, y::ForwardDiff.Dual) = x^y

sse(x::Number) = x^2
sse(x::ForwardDiff.Dual) = sse(ForwardDiff.value(x)) + sum(sse, ForwardDiff.partials(x))
totallength(x::Number) = 1
function totallength(x::ForwardDiff.Dual)
    totallength(ForwardDiff.value(x)) + sum(totallength, ForwardDiff.partials(x))
end
totallength(x::AbstractArray) = sum(totallength, x)

@inline ODE_DEFAULT_NORM(u::ForwardDiff.Dual, ::Any) = sqrt(sse(u))
@inline function ODE_DEFAULT_NORM(u::AbstractArray{<:ForwardDiff.Dual}, t::Any)
    sqrt(sum(sse, u) / totallength(u))
end
@inline ODE_DEFAULT_NORM(u::ForwardDiff.Dual, ::ForwardDiff.Dual) = sqrt(sse(u))
@inline function ODE_DEFAULT_NORM(u::AbstractArray{<:ForwardDiff.Dual}, ::ForwardDiff.Dual)
    sqrt(sum(sse, u) / totallength(u))
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
