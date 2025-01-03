module DiffEqBaseGTPSAExt

if isdefined(Base, :get_extension)
    using DiffEqBase
    import DiffEqBase: value
    using GTPSA
else
    using ..DiffEqBase
    import ..DiffEqBase: value
    using ..GTPSA
end

value(x::TPS) = scalar(x);
value(::Type{TPS{T}}) where {T} = T


end