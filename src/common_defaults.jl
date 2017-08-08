@inline UNITLESS_ABS2(x) = Real(abs2(x)/(typeof(x)(one(x))*typeof(x)(one(x))))
@inline ODE_DEFAULT_NORM(u::AbstractArray) = sqrt(sum(UNITLESS_ABS2,u) / length(u))
@inline ODE_DEFAULT_NORM(u::AbstractArray{T,N}) where {T<:AbstractArray,N} = sqrt(sum(ODE_DEFAULT_NORM,u) / length(u))
@inline ODE_DEFAULT_NORM(u) = norm(u)
@inline ODE_DEFAULT_ISOUTOFDOMAIN(t,u) = false
@inline ODE_DEFAULT_PROG_MESSAGE(dt,t,u) = "dt="*string(dt)*"\nt="*string(t)*"\nmax u="*string(maximum(abs.(u)))
@inline ODE_DEFAULT_UNSTABLE_CHECK(dt,t,u) = false
(p::typeof(ODE_DEFAULT_UNSTABLE_CHECK))(dt,t,u::AbstractFloat) = isnan(u)
(p::typeof(ODE_DEFAULT_UNSTABLE_CHECK))(dt,t,u::AbstractArray{T}) where {T<:AbstractFloat} = any(isnan,u)
(p::typeof(ODE_DEFAULT_UNSTABLE_CHECK))(dt,t,u::ArrayPartition) =
                                                 any(any(isnan,x) for x in u.x)
