type SDEProblem{uType,tType,isinplace,isinplaceNoise,NoiseClass,F,F2,F3,C,ND} <: AbstractSDEProblem{uType,tType,isinplace,NoiseClass,ND}
  f::F
  g::F2
  u0::uType
  tspan::Tuple{tType,tType}
  noise::NoiseProcess{NoiseClass,isinplaceNoise,F3}
  callback::C
  noise_dim::ND
end

function SDEProblem{NoiseClass,inplace,F3}(f,g,u0,tspan;iip = isinplace(f,3),
                    noise::NoiseProcess{NoiseClass,inplace,F3}= iip ? INPLACE_WHITE_NOISE : WHITE_NOISE,
                    callback = CallbackSet(),noise_dim = nothing)
  SDEProblem{typeof(u0),eltype(tspan),iip,inplace,NoiseClass,typeof(f),typeof(g),F3,typeof(callback),typeof(noise_dim)}(f,g,u0,tspan,noise,callback,noise_dim)
end
