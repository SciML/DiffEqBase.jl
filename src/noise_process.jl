type NoiseProcess{class,inplace,F}
  noise_func::F
end

(n::NoiseProcess)(a) = n.noise_func(a)
(n::NoiseProcess)(a...) = n.noise_func(a...)
(n::NoiseProcess)(a,b) = n.noise_func(a,b)

const WHITE_NOISE = NoiseProcess{:Diagonal,false,typeof(randn)}(randn)
const INPLACE_WHITE_NOISE = NoiseProcess{:Diagonal,true,typeof(randn!)}(randn!)

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
  NoiseProcess{:Commutative,true,typeof(noise_func!)}(noise_func!)
end

isinplace{class,inplace,F}(n::NoiseProcess{class,inplace,F}) = inplace
