module DiffEqBaseGTPSAExt

if isdefined(Base, :get_extension)
    using DiffEqBase
    import DiffEqBase: value, ODE_DEFAULT_NORM
    using GTPSA
else
    using ..DiffEqBase
    import ..DiffEqBase: value, ODE_DEFAULT_NORM
    using ..GTPSA
end

value(x::TPS) = scalar(x)
value(::Type{TPS{T}}) where {T} = T

ODE_DEFAULT_NORM(u::TPS, t) = normTPS(u)
ODE_DEFAULT_NORM(f::F, u::TPS, t) where {F} = normTPS(f(u))

function ODE_DEFAULT_NORM(u::AbstractArray{TPS{T}}, t) where {T}
    x = zero(real(T))
    @inbounds @fastmath for ui in u
        x += normTPS(ui)^2
    end
    Base.FastMath.sqrt_fast(x / max(length(u), 1))
end

function ODE_DEFAULT_NORM(f::F, u::AbstractArray{TPS{T}}, t) where {F, T}
    x = zero(real(T))
    @inbounds @fastmath for ui in u
        x += normTPS(f(ui))^2
    end
    Base.FastMath.sqrt_fast(x / max(length(u), 1))
end

end