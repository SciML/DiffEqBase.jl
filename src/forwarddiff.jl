promote_dual(::Type{T},::Type{T2}) where {T,T2} = T
promote_dual(::Type{T},::Type{T2}) where {T<:ForwardDiff.Dual,T2} = T
promote_dual(::Type{T},::Type{T2}) where {T<:ForwardDiff.Dual,T2<:ForwardDiff.Dual} = T
promote_dual(::Type{T},::Type{T2}) where {T,T2<:ForwardDiff.Dual} = T2

anyeltypedual(x) = mapreduce(y->anyeltypedual(getproperty(x,y)),promote_dual,propertynames(x))
anyeltypedual(x::Union{String,Symbol}) = Any
anyeltypedual(::Type{T}) where T = T
anyeltypedual(x::Number) = typeof(x)
anyeltypedual(x::Union{AbstractArray{T},Set{T}}) where T<:Number = T
anyeltypedual(x::Union{AbstractArray,Set}) = mapreduce(anyeltypedual,promote_dual,x)
anyeltypedual(x::Tuple) = mapreduce(anyeltypedual,promote_dual,x)
anyeltypedual(x::Union{Dict,NamedTuple}) = mapreduce(anyeltypedual,promote_dual,values(x))

function promote_u0(u0,p,t0) 
  T = anyeltypedual(p)
  if T <: ForwardDiff.Dual && !(eltype(u0) <: ForwardDiff.Dual)
    T.(u0)
  else
    u0
  end
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
