@inline UNITLESS_ABS2(x) = Real(abs2(x)/(typeof(x)(one(x))*typeof(x)(one(x))))
@inline ODE_DEFAULT_NORM(u::Union{AbstractFloat,Complex}) = abs(u)
@inline ODE_DEFAULT_NORM(u::Array{T}) where T<:Union{AbstractFloat,Complex} =
                                         @fastmath sqrt(sum(abs2,u) / length(u))
@inline ODE_DEFAULT_NORM(u::AbstractArray) = sqrt(sum(UNITLESS_ABS2,u) / length(u))
@inline ODE_DEFAULT_NORM(u::AbstractArray{T,N}) where {T<:AbstractArray,N} = sqrt(sum(ODE_DEFAULT_NORM,u) / length(u))
@inline ODE_DEFAULT_NORM(u) = norm(u)
@inline ODE_DEFAULT_ISOUTOFDOMAIN(u,p,t) = false
@inline ODE_DEFAULT_PROG_MESSAGE(dt,u,p,t) =
           "dt="*string(dt)*"\nt="*string(t)*"\nmax u="*string(maximum(abs.(u)))
@inline ODE_DEFAULT_UNSTABLE_CHECK(dt,u,p,t) = false
@inline ODE_DEFAULT_UNSTABLE_CHECK(dt,u,p,t::AbstractFloat) = isnan(u)
@inline ODE_DEFAULT_UNSTABLE_CHECK(dt,u,p,t::AbstractArray{T}) where
                                               {T<:AbstractFloat} = any(isnan,u)
@inline ODE_DEFAULT_UNSTABLE_CHECK(dt,u,p,t::ArrayPartition) =
                                                 any(any(isnan,x) for x in u.x)
