"""
construct_correlated_noisefunc(Γ::AbstractArray)

Takes in a constant Covariance matrix Γ and spits out the noisefunc.
"""
function construct_correlated_noisefunc(Γ::AbstractArray)
  γ = svdfact(Γ)
  A = γ[:U]*Diagonal(√γ[:S])
  b = Vector{eltype(Γ)}(size(Γ,1))
  noise_func! = function (a,integrator)
    randn!(b)
    A_mul_B!(a,A,b)
  end
  NoiseProcess{:White,true,typeof(noise_func!)}(noise_func!)
end
