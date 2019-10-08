@inline UNITLESS_ABS2(x) = real(abs2(x))
@inline ODE_DEFAULT_NORM(u::Union{AbstractFloat,Complex},t) = abs(u)
@inline ODE_DEFAULT_NORM(u::Array{T},t) where T<:Union{AbstractFloat,Complex} =
                                         sqrt(real(sum(abs2,u)) / length(u))
@inline ODE_DEFAULT_NORM(u::StaticArray{T},t) where T<:Union{AbstractFloat,Complex} =
                                            sqrt(real(sum(abs2,u)) / length(u))
@inline ODE_DEFAULT_NORM(u::RecursiveArrayTools.AbstractVectorOfArray,t) =
                                        sum(sqrt(real(sum(UNITLESS_ABS2,_u)) / length(_u)) for _u in u.u)
@inline ODE_DEFAULT_NORM(u::AbstractArray,t) = sqrt(real(sum(UNITLESS_ABS2,u)) / length(u))
@inline ODE_DEFAULT_NORM(u::AbstractArray{T,N},t) where {T<:AbstractArray,N} = sqrt(real(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip(u,t))) / length(u))
@inline ODE_DEFAULT_NORM(u,t) = norm(u)

@inline ODE_DEFAULT_ISOUTOFDOMAIN(u,p,t) = false
@inline ODE_DEFAULT_PROG_MESSAGE(dt,u,p,t) =
           "dt="*string(dt)*"\nt="*string(t)*"\nmax u="*string(maximum(abs.(u)))
@inline ODE_DEFAULT_UNSTABLE_CHECK(dt,u,p,t) = false
@inline ODE_DEFAULT_UNSTABLE_CHECK(dt,u::AbstractFloat,p,t) = isnan(u)
@inline ODE_DEFAULT_UNSTABLE_CHECK(dt,u::Float64,p,t) =
                                                any(x->(isnan(x) || x>1e50),u)
@inline ODE_DEFAULT_UNSTABLE_CHECK(dt,u::AbstractArray{T},p,t) where
                                    {T<:AbstractFloat} = any(isnan,u)
@inline ODE_DEFAULT_UNSTABLE_CHECK(dt,u::ArrayPartition,p,t) =
                                                 any(any(isnan,x) for x in u.x)
