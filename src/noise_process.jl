type NoiseProcess{class,F}
  noise_func::F
end

const WHITE_NOISE = NoiseProcess{:Diagonal,typeof(randn)}(randn)

"""
construct_correlated_noisefunc(Γ::AbstractArray)

Takes in a constant Covariance matrix Γ and spits out the noisefunc.
"""
function construct_correlated_noisefunc(Γ::AbstractArray)
  γ = svdfact(Γ)
  A = γ[:U]*Diagonal(√γ[:S])
  noise_func = (N...) -> A*randn(N...)
  NoiseProcess{:Commutative,typeof(noise_func)}(noise_func)
end
