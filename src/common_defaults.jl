abs2_and_sum(x,y) = reduce(Base.add_sum,x,init=zero(real(value(eltype(x))))) +
                    reduce(Base.add_sum,y,init=zero(real(value(eltype(y)))))
@inline UNITLESS_ABS2(x::Number) = abs2(x)
@inline UNITLESS_ABS2(x::AbstractArray) = mapreduce(UNITLESS_ABS2,abs2_and_sum, x, init=zero(real(value(eltype(x)))))
@inline UNITLESS_ABS2(x::RecursiveArrayTools.ArrayPartition) = mapreduce(UNITLESS_ABS2, abs2_and_sum, x.x, init=zero(real(value(eltype(x)))))

@inline recursive_length(u::AbstractArray{<:Number}) = length(u)
@inline recursive_length(u::Number) = length(u)
@inline recursive_length(u::AbstractArray{<:AbstractArray}) = sum(recursive_length, u)
@inline recursive_length(u::RecursiveArrayTools.ArrayPartition) = sum(recursive_length, u.x)
@inline recursive_length(u::RecursiveArrayTools.VectorOfArray) = sum(recursive_length, u.u)
@inline function recursive_length(u::AbstractArray{<:StaticArray{Size, <:Number}}) where {Size}
  prod(Size(eltype(u))) * length(u)
end

@inline ODE_DEFAULT_NORM(u::Union{AbstractFloat,Complex},t) = @fastmath abs(u)

@inline function ODE_DEFAULT_NORM(u::Array{T},t) where T<:Union{AbstractFloat,Complex}
    x = abs2(u[1])
    @inbounds for i in 2:length(u)
        x += abs2(u[i])
    end
    sqrt(real(x) / length(u))
end

@inline ODE_DEFAULT_NORM(u::StaticArrays.StaticArray{T},t) where T<:Union{AbstractFloat,Complex} =
                                            sqrt(real(sum(abs2,u)) / length(u))

@inline ODE_DEFAULT_NORM(u::AbstractArray,t) = sqrt(UNITLESS_ABS2(u)/recursive_length(u))
@inline ODE_DEFAULT_NORM(u,t) = norm(u)

@inline ODE_DEFAULT_ISOUTOFDOMAIN(u,p,t) = false
@inline function ODE_DEFAULT_PROG_MESSAGE(dt,u::Array,p,t)
    tmp = u[1]
    for i in eachindex(u)
        tmp = ifelse(abs(u[i]) > abs(tmp),u[i],tmp)
    end
    "dt="*string(dt)*"\nt="*string(t)*"\nmax u="*string(tmp)
end
@inline ODE_DEFAULT_PROG_MESSAGE(dt,u,p,t) =
           "dt="*string(dt)*"\nt="*string(t)*"\nmax u="*string(maximum(abs.(u)))

@inline NAN_CHECK(x::Number) = isnan(x)
@inline NAN_CHECK(x::Float64) = isnan(x) || (x>1e50)
@inline NAN_CHECK(x::Enum) = false
@inline NAN_CHECK(x::AbstractArray) = any(NAN_CHECK, x)
@inline NAN_CHECK(x::RecursiveArrayTools.ArrayPartition) = any(NAN_CHECK, x.x)
@inline NAN_CHECK(x::SparseArrays.AbstractSparseMatrixCSC) = any(NAN_CHECK, nonzeros(x))
    
@inline ODE_DEFAULT_UNSTABLE_CHECK(dt,u,p,t) = false
@inline ODE_DEFAULT_UNSTABLE_CHECK(dt,u::Union{Number,AbstractArray},p,t) = NAN_CHECK(u)
