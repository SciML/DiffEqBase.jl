module DiffEqBaseUnitfulExt

if isdefined(Base, :get_extension)
    using DiffEqBase
    import DiffEqBase: value
    using Unitful
else
    using ..DiffEqBase
    import ..DiffEqBase: value
    using ..Unitful
end
# Support adaptive errors should be errorless for exponentiation
value(x::Type{Unitful.AbstractQuantity{T, D, U}}) where {T, D, U} = T
value(x::Unitful.AbstractQuantity) = x.val
@inline function DiffEqBase.ODE_DEFAULT_NORM(u::AbstractArray{
        <:Unitful.AbstractQuantity,
        N,
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

_rate_prototype(u, t, onet) = u / unit(t)
end
