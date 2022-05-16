promote_u0(u0::AbstractArray{<:ForwardDiff.Dual},p::AbstractArray{<:ForwardDiff.Dual},t0) = u0
promote_u0(u0,p::AbstractArray{<:ForwardDiff.Dual},t0) = eltype(p).(u0)
promote_u0(u0,p::NTuple{N,<:ForwardDiff.Dual},t0) where N = eltype(p).(u0)
promote_u0(u0,p::ForwardDiff.Dual,t0) where N = eltype(p).(u0)

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
