@inline wiener_randn() = randn()
@inline wiener_randn(x...) = randn(x...)
@inline wiener_randn!(x...) = randn!(x...)
@inline wiener_randn{T<:Number}(::Type{Complex{T}}) = (randn(T)+im*randn(T))/sqrt(2)
@inline wiener_randn{T<:Number}(::Type{Complex{T}},x...) = (randn(T,x...)+im*randn(T,x...))/sqrt(2)
@inline wiener_randn{T<:Number}(y::AbstractRNG,::Type{Complex{T}},x...) = (randn(y,T,x...)+im*randn(y,T,x...))/sqrt(2)
@inline wiener_randn{T<:Number}(y::AbstractRNG,::Type{Complex{T}}) = (randn(y,T)+im*randn(y,T))/sqrt(2)
@inline function wiener_randn!{T<:Number}(y::AbstractRNG,x::AbstractArray{Complex{T}})
  for i in eachindex(x)
    x[i] = (randn(y,T)+im*randn(y,T))/sqrt(2)
  end
end
@inline function wiener_randn!{T<:Number}(x::AbstractArray{Complex{T}})
  for i in eachindex(x)
    x[i] = (randn(T)+im*randn(T))/sqrt(2)
  end
end

type NoiseProcess{class,inplace,F}
  noise_func::F
end

(n::NoiseProcess)(a) = n.noise_func(a)
(n::NoiseProcess)(a...) = n.noise_func(a...)
(n::NoiseProcess)(a,b) = n.noise_func(a,b)

const WHITE_NOISE = NoiseProcess{:White,false,typeof(wiener_randn)}(wiener_randn)
const INPLACE_WHITE_NOISE = NoiseProcess{:White,true,typeof(wiener_randn!)}(wiener_randn!)

"""
construct_correlated_noisefunc(Γ::AbstractArray)

Takes in a constant Covariance matrix Γ and spits out the noisefunc.
"""
function construct_correlated_noisefunc(Γ::AbstractArray)
  γ = svdfact(Γ)
  A = γ[:U]*Diagonal(√γ[:S])
  b = Vector{eltype(Γ)}(size(Γ,1))
  noise_func! = function (a)
    randn!(b)
    A_mul_B!(a,A,b)
  end
  NoiseProcess{:White,true,typeof(noise_func!)}(noise_func!)
end

isinplace{class,inplace,F}(n::NoiseProcess{class,inplace,F}) = inplace
noise_class{class,inplace,F}(n::NoiseProcess{class,inplace,F}) = class
