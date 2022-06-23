"""
    promote_dual(::Type{T},::Type{T2})


Is like the number promotion system, but always prefers a dual number type above
anything else. For higher order differentiation, it returns the most dualiest of
them all. This is then used to promote `u0` into the suspected highest differentiation
space for solving the equation.
"""
promote_dual(::Type{T},::Type{T2}) where {T,T2} = T
promote_dual(::Type{T},::Type{T2}) where {T<:ForwardDiff.Dual,T2} = T
promote_dual(::Type{T},::Type{T2}) where {T<:ForwardDiff.Dual,T2<:ForwardDiff.Dual} = T
promote_dual(::Type{T},::Type{T2}) where {T,T2<:ForwardDiff.Dual} = T2

promote_dual(::Type{T},::Type{T2}) where {T3,T4,V,V2<:ForwardDiff.Dual,N,N2,
             T<:ForwardDiff.Dual{T3,V,N},
             T2<:ForwardDiff.Dual{T4,V2,N2}} = T2
promote_dual(::Type{T},::Type{T2}) where {T3,T4,V<:ForwardDiff.Dual,V2,N,N2,
             T<:ForwardDiff.Dual{T3,V,N},
             T2<:ForwardDiff.Dual{T4,V2,N2}} = T
promote_dual(::Type{T},::Type{T2}) where {
             T3,V<:ForwardDiff.Dual,V2<:ForwardDiff.Dual,N,
             T<:ForwardDiff.Dual{T3,V,N},
             T2<:ForwardDiff.Dual{T3,V2,N}} = ForwardDiff.Dual{T3,promote_dual(V,V2),N}

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
function anyeltypedual(x) 
  if propertynames(x) === ()
    Any
  else
    mapreduce(y-> !isdefined(x,y) ? Any : anyeltypedual(getproperty(x,y)),promote_dual,propertynames(x))
  end
end

Base.@pure anyeltypedual(::Type{T}) where T = hasproperty(T,:parameters) ? mapreduce(anyeltypedual,promote_dual,T.parameters;init=Any) : T
anyeltypedual(::Type{T}) where T<:Union{ForwardDiff.Dual} = T
anyeltypedual(::Type{T}) where T<:Union{AbstractArray,Set} = anyeltypedual(eltype(T))
Base.@pure function anyeltypedual(::Type{T}) where T<:Union{NTuple}
  if isconcretetype(eltype(T)) 
    return eltype(T) 
  end
  if isempty(T.parameters)
    Any
  else
    mapreduce(anyeltypedual,promote_dual,T.parameters;init=Any)
  end
end

# Any in this context just means not Dual
anyeltypedual(x::SciMLBase.NullParameters) = Any
anyeltypedual(x::Number) = anyeltypedual(typeof(x))
anyeltypedual(x::Union{String,Symbol}) = typeof(x)
anyeltypedual(x::Union{Array{T},AbstractArray{T},Set{T}}) where T<:Union{Number,Symbol,String} = anyeltypedual(T)
anyeltypedual(x::Union{Array{T},AbstractArray{T},Set{T}}) where {N,T<:Union{AbstractArray{<:Number},Set{<:Number},NTuple{N,<:Number}}} = anyeltypedual(eltype(x))

# Try to avoid this dispatch because it can lead to type inference issues when !isconcrete(eltype(x))
function anyeltypedual(x::AbstractArray) 
  if isconcretetype(eltype(x))
    anyeltypedual(eltype(x))
  elseif !isempty(x) && all(i->isassigned(x,i),1:length(x))
    mapreduce(anyeltypedual,promote_dual,x)
  else
    # This fallback to Any is required since otherwise we cannot handle `undef` in all cases
    #  misses cases of 
    Any
  end
end

function anyeltypedual(x::Set) 
  if isconcretetype(eltype(x))
    anyeltypedual(eltype(x))
  else
    # This fallback to Any is required since otherwise we cannot handle `undef` in all cases
    Any
  end
end

function anyeltypedual(x::Tuple)
  # Handle the empty tuple case separately for inference and to avoid mapreduce error
  if x === ()
    Any
  else
    mapreduce(anyeltypedual,promote_dual,x)
  end
end
anyeltypedual(x::Dict) = isempty(x) ? eltype(values(x)) : mapreduce(anyeltypedual,promote_dual,values(x))
anyeltypedual(x::NamedTuple) = isempty(x) ? Any : mapreduce(anyeltypedual,promote_dual,values(x))

function promote_u0(u0,p,t0) 
  if !(eltype(u0) <: ForwardDiff.Dual)
    T = anyeltypedual(p)
    if T <: ForwardDiff.Dual 
      return T.(u0)
    end
  end
  u0
end

function promote_tspan(u0::AbstractArray{<:ForwardDiff.Dual},p,tspan::Tuple{<:ForwardDiff.Dual,<:ForwardDiff.Dual},prob,kwargs)
  return tspan
end

function promote_tspan(u0::AbstractArray{<:ForwardDiff.Dual},p,tspan,prob,kwargs)
  if (haskey(kwargs,:callback) && has_continuous_callback(kwargs[:callback])) ||
     (haskey(prob.kwargs,:callback) && has_continuous_callback(prob.kwargs[:callback]))

    return eltype(u0).(tspan)
  else
    return tspan
  end
end

value(x::Type{ForwardDiff.Dual{T,V,N}}) where {T,V,N} = V
value(x::ForwardDiff.Dual) = value(ForwardDiff.value(x))

@inline fastpow(x::ForwardDiff.Dual, y::ForwardDiff.Dual) = x^y

sse(x::Number) = x^2
sse(x::ForwardDiff.Dual) = sse(ForwardDiff.value(x)) + sum(sse, ForwardDiff.partials(x))
totallength(x::Number) = 1
totallength(x::ForwardDiff.Dual) = totallength(ForwardDiff.value(x)) + sum(totallength, ForwardDiff.partials(x))
totallength(x::AbstractArray) = sum(totallength,x)

@inline ODE_DEFAULT_NORM(u::ForwardDiff.Dual,::Any) = sqrt(sse(u))
@inline ODE_DEFAULT_NORM(u::AbstractArray{<:ForwardDiff.Dual},t::Any) = sqrt(sum(sse,u) / totallength(u))
@inline ODE_DEFAULT_NORM(u::ForwardDiff.Dual,::ForwardDiff.Dual) = sqrt(sse(u))
@inline ODE_DEFAULT_NORM(u::AbstractArray{<:ForwardDiff.Dual},::ForwardDiff.Dual) = sqrt(sum(sse,u) / totallength(u))

if !hasmethod(nextfloat, Tuple{ForwardDiff.Dual})
  # Type piracy. Should upstream
  Base.nextfloat(d::ForwardDiff.Dual{T,V,N}) where {T,V,N} = ForwardDiff.Dual{T}(nextfloat(d.value), d.partials)
  Base.prevfloat(d::ForwardDiff.Dual{T,V,N}) where {T,V,N} = ForwardDiff.Dual{T}(prevfloat(d.value), d.partials)
end

# bisection(f, tup::Tuple{T,T}, t_forward::Bool) where {T<:ForwardDiff.Dual} = find_zero(f, tup, Roots.AlefeldPotraShi())
