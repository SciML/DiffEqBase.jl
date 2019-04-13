import Base.Broadcast: _broadcast_getindex, preprocess, preprocess_args, Broadcasted, broadcast_unalias, combine_axes, broadcast_axes, broadcast_shape, check_broadcast_axes, check_broadcast_shape, throwdm
import Base: copyto!, tail, axes
struct DiffEqBC{T}
    x::T
end
@inline axes(b::DiffEqBC) = axes(b.x)
Base.@propagate_inbounds _broadcast_getindex(b::DiffEqBC, i) = _broadcast_getindex(b.x, i)
Base.@propagate_inbounds _broadcast_getindex(b::DiffEqBC{<:AbstractArray{<:Any,0}}, i) = b.x[]
Base.@propagate_inbounds _broadcast_getindex(b::DiffEqBC{<:AbstractVector}, i) = b.x[i[1]]
Base.@propagate_inbounds _broadcast_getindex(b::DiffEqBC{<:AbstractArray}, i) = b.x[i]
diffeqbc(x::Array) = DiffEqBC(x)
diffeqbc(x) = x

# Ensure inlining
@inline combine_axes(A, B) = broadcast_shape(broadcast_axes(A), broadcast_axes(B)) # Julia 1.0 compatible
@inline check_broadcast_axes(shp, A::Union{Number, Array, Broadcasted}) = check_broadcast_shape(shp, axes(A))

@inline preprocess(f, dest, bc::Broadcasted{Style}) where {Style} = Broadcasted{Style}(bc.f, preprocess_args(f, dest, bc.args), bc.axes)
preprocess(f, dest, x) = f(broadcast_unalias(dest, x))

@inline preprocess_args(f, dest, args::Tuple) = (preprocess(f, dest, args[1]), preprocess_args(f, dest, tail(args))...)
@inline preprocess_args(f, dest, args::Tuple{Any}) = (preprocess(f, dest, args[1]),)
preprocess_args(f, dest, args::Tuple{}) = ()

@inline function copyto!(dest::DiffEqBC, bc::Broadcasted)
    axes(dest) == axes(bc) || throwdm(axes(dest), axes(bc))
    bcs′ = preprocess(diffeqbc, dest, bc)
    dest′ = dest.x
    @simd ivdep for I in eachindex(bcs′)
        @inbounds dest′[I] = bcs′[I]
    end
    return dest′ # return the original array without the wrapper
end

import Base.Broadcast: broadcasted, broadcastable, combine_styles
map_nostop(f, t::Tuple{})              = ()
map_nostop(f, t::Tuple{Any,})          = (f(t[1]),)
map_nostop(f, t::Tuple{Any, Any})      = (f(t[1]), f(t[2]))
map_nostop(f, t::Tuple{Any, Any, Any}) = (f(t[1]), f(t[2]), f(t[3]))
map_nostop(f, t::Tuple)                = (Base.@_inline_meta; (f(t[1]), map_nostop(f,tail(t))...))
@inline function broadcasted(f::Union{typeof(*), typeof(+), typeof(muladd)}, arg1, arg2, args...)
    arg1′ = broadcastable(arg1)
    arg2′ = broadcastable(arg2)
    args′ = map_nostop(broadcastable, args)
    broadcasted(combine_styles(arg1′, arg2′, args′...), f, arg1′, arg2′, args′...)
end

macro ..(x)
    expr = Base.Broadcast.__dot__(x)
    if expr.head == :(.=)
      dest = expr.args[1]
      expr.args[1] = :(DiffEqBase.diffeqbc($dest))
    end
    return esc(expr)
end
