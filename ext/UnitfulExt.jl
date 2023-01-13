module UnitfulExt

using DiffEqBase
isdefined(Base, :get_extension) ? (using Unitful) : (using ..Unitful)

# Support adaptive errors should be errorless for exponentiation
DiffEqBase.value(x::Type{Unitful.AbstractQuantity{T, D, U}}) where {T, D, U} = T
DiffEqBase.value(x::Unitful.AbstractQuantity) = x.val
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::AbstractArray{<:Unitful.AbstractQuantity, N
                                                              },
                                             t) where {N}
    sqrt(sum(x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]),
             zip((value(x) for x in u), Iterators.repeated(t))) / length(u))
end
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::Array{<:Unitful.AbstractQuantity, N},
                                             t) where {N}
    sqrt(sum(x -> DiffEqBase.ODE_DEFAULT_NORM(x[1], x[2]),
             zip((value(x) for x in u), Iterators.repeated(t))) / length(u))
end
@inline DiffEqBase.ODE_DEFAULT_NORM(u::Unitful.AbstractQuantity, t) = abs(value(u))
@inline function DiffEqBase.UNITLESS_ABS2(x::Unitful.AbstractQuantity)
    real(abs2(x) / oneunit(x) * oneunit(x))
end

end
