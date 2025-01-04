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

ODE_DEFAULT_NORM(u::TPS, t) = @fastmath abs(value(u))
ODE_DEFAULT_NORM(f::F, u::TPS, t) where {F} = @fastmath abs(f(value(u)))

function ODE_DEFAULT_NORM(u::AbstractArray{TPS{T}}, t) where {T <: Union{AbstractFloat, Complex}}
    x = zero(real(T))
    @inbounds @fastmath for ui in u
        x += abs2(value(ui))
    end
    Base.FastMath.sqrt_fast(x / max(length(u), 1))
end

end