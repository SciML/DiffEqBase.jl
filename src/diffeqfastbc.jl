import Base.Broadcast: _broadcast_getindex, preprocess, preprocess_args, Broadcasted, broadcast_unalias, combine_axes, broadcast_shape, check_broadcast_axes, check_broadcast_shape
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

@inline combine_axes(A, B, C...) = broadcast_shape(axes(A), combine_axes(B, C...))
@inline combine_axes(A, B) = broadcast_shape(axes(A), axes(B))
@inline check_broadcast_axes(shp, A) = check_broadcast_shape(shp, axes(A))

@inline preprocess(f, dest, bc::Broadcasted{Style}) where {Style} = Broadcasted{Style}(bc.f, preprocess_args(f, dest, bc.args), bc.axes)
preprocess(f, dest, x) = f(broadcast_unalias(dest, x))

@inline preprocess_args(f, dest, args::Tuple) = (preprocess(f, dest, args[1]), preprocess_args(f, dest, tail(args))...)
@inline preprocess_args(f, dest, args::Tuple{Any}) = (preprocess(f, dest, args[1]),)
preprocess_args(f, dest, args::Tuple{}) = ()

@static if VERSION >= v"1.2.0"
@eval Base.getindex(A::DiffEqBC, i1::Int) =
    (Base.@_inline_meta; Core.const_arrayref($(Expr(:boundscheck)), A.x, i1))
@eval Base.getindex(A::DiffEqBC, i1::Int, i2::Int, I::Int...) =
  (Base.@_inline_meta; Core.const_arrayref($(Expr(:boundscheck)), A.x, i1, i2, I...))
macro aliasscope(body)
    sym = gensym()
    quote
        $(Expr(:aliasscope))
        $sym = $(esc(body))
        $(Expr(:popaliasscope))
        $sym
    end
end
end

@inline function copyto!(dest::DiffEqBC, bc::Broadcasted)
    axes(dest) == axes(bc) || throwdm(axes(dest), axes(bc))
    bcs′ = preprocess(diffeqbc, dest, bc)
    dest′ = dest.x
    @static if VERSION >= v"1.2.0"
        @aliasscope @simd for I in eachindex(bcs′)
            @inbounds dest′[I] = bcs′[I]
        end
    else
        @simd for I in eachindex(bcs′)
            @inbounds dest′[I] = bcs′[I]
        end
    end
    return dest
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
