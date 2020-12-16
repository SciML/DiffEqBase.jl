@inline UNITLESS_ABS2(x::Number) = abs2(x)
@inline UNITLESS_ABS2(x::AbstractArray) = sum(UNITLESS_ABS2, x)
@inline UNITLESS_ABS2(x::ArrayPartition) = sum(UNITLESS_ABS2, x.x)

@inline recursive_length(u::AbstractArray{<:Number}) = length(u)
@inline recursive_length(u::Number) = 1
@inline recursive_length(u::AbstractArray{<:AbstractArray}) = sum(recursive_length, u)
@inline recursive_length(u::ArrayPartition) = sum(recursive_length, u.x)
@inline recursive_length(u::VectorOfArray) = sum(recursive_length, u.u)

@inline ODE_DEFAULT_NORM(u::Union{AbstractFloat,Complex},t) = @fastmath abs(u)
@inline ODE_DEFAULT_NORM(u::Array{T},t) where T<:Union{AbstractFloat,Complex} =
                                         sqrt(real(sum(abs2,u)) / length(u))
@inline ODE_DEFAULT_NORM(u::StaticArray{T},t) where T<:Union{AbstractFloat,Complex} =
                                            sqrt(real(sum(abs2,u)) / length(u))

@inline ODE_DEFAULT_NORM(u::AbstractArray,t) = sqrt(UNITLESS_ABS2(u)/recursive_length(u))
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
