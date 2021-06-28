@inline UNITLESS_ABS2(x::Number) = abs2(x)
@inline UNITLESS_ABS2(x::AbstractArray) = sum(UNITLESS_ABS2, x)
@inline UNITLESS_ABS2(x::RecursiveArrayTools.ArrayPartition) = sum(UNITLESS_ABS2, x.x)
Base.mapreduce_empty(::typeof(UNITLESS_ABS2), op, T) = abs2(Base.reduce_empty(op, T))

@inline recursive_length(u::AbstractArray{<:Number}) = length(u)
@inline recursive_length(u::Number) = length(u)
@inline recursive_length(u::AbstractArray{<:AbstractArray}) = sum(recursive_length, u)
@inline recursive_length(u::RecursiveArrayTools.ArrayPartition) = sum(recursive_length, u.x)
@inline recursive_length(u::RecursiveArrayTools.VectorOfArray) = sum(recursive_length, u.u)

@inline ODE_DEFAULT_NORM(u::Union{AbstractFloat,Complex},t) = @fastmath abs(u)
@inline ODE_DEFAULT_NORM(u::Array{T},t) where T<:Union{AbstractFloat,Complex} =
                                         sqrt(real(sum(abs2,u)) / length(u))
@inline ODE_DEFAULT_NORM(u::StaticArrays.StaticArray{T},t) where T<:Union{AbstractFloat,Complex} =
                                            sqrt(real(sum(abs2,u)) / length(u))

@inline ODE_DEFAULT_NORM(u::AbstractArray,t) = sqrt(UNITLESS_ABS2(u)/recursive_length(u))
@inline ODE_DEFAULT_NORM(u,t) = norm(u)

@inline ODE_DEFAULT_ISOUTOFDOMAIN(u,p,t) = false
@inline ODE_DEFAULT_PROG_MESSAGE(dt,u,p,t) =
           "dt="*string(dt)*"\nt="*string(t)*"\nmax u="*string(maximum(abs.(u)))

@inline NAN_CHECK(x::Number) = isnan(x)
@inline NAN_CHECK(x::Float64) = isnan(x) || (x>1e50)
@inline NAN_CHECK(x::Enum) = false
@inline NAN_CHECK(x::AbstractArray) = any(NAN_CHECK, x)
@inline NAN_CHECK(x::RecursiveArrayTools.ArrayPartition) = any(NAN_CHECK, x.x)
@inline NAN_CHECK(x::DEDataArray) = NAN_CHECK(x.x)


@inline ODE_DEFAULT_UNSTABLE_CHECK(dt,u,p,t) = false
@inline ODE_DEFAULT_UNSTABLE_CHECK(dt,u::Union{Number,AbstractArray},p,t) = NAN_CHECK(u)
