type RODEProblem{uType,tType,isinplace,NoiseClass,isinplaceNoise,F,C,MM,F3,ND} <: AbstractRODEProblem{uType,tType,isinplace,NoiseClass,ND}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  noise::NoiseProcess{NoiseClass,isinplaceNoise,F3}
  callback::C
  mass_matrix::MM
  rand_prototype::ND
end

function RODEProblem{NoiseClass,inplace,F3}(f,u0,tspan; iip = isinplace(f,4),
                     noise::NoiseProcess{NoiseClass,inplace,F3}= iip ? INPLACE_WHITE_NOISE : WHITE_NOISE,
                     callback=nothing,rand_prototype = nothing,mass_matrix=I)
  RODEProblem{typeof(u0),eltype(tspan),iip,NoiseClass,inplace,typeof(f),typeof(callback),typeof(mass_matrix),F3,typeof(rand_prototype)}(f,u0,tspan,noise,callback,mass_matrix,rand_prototype)
end
