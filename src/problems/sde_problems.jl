type SDEProblem{uType,tType,isinplace,isinplaceNoise,NoiseClass,F,F2,F3} <: AbstractSDEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3}
  f::F
  g::F2
  u0::uType
  tspan::Tuple{tType,tType}
  noise::NoiseProcess{NoiseClass,isinplaceNoise,F3}
end

type SDETestProblem{uType,tType,isinplace,isinplaceNoise,NoiseClass,F,F2,F3} <: AbstractSDETestProblem{uType,tType,isinplace,NoiseClass,F,F2,F3}
  f::F
  g::F2
  u0::uType
  analytic::Base.Callable
  tspan::Tuple{tType,tType}
  noise::NoiseProcess{NoiseClass,isinplaceNoise,F3}
end

function SDEProblem{NoiseClass,inplace,F3}(f,g,u0,tspan;iip = isinplace(f,3),
                    noise::NoiseProcess{NoiseClass,inplace,F3}= iip ? INPLACE_WHITE_NOISE : WHITE_NOISE)
  SDEProblem{typeof(u0),eltype(tspan),iip,inplace,NoiseClass,typeof(f),typeof(g),F3}(f,g,u0,tspan,noise)
end

function SDETestProblem{NoiseClass,inplace,F3}(f,g,u0,analytic,tspan=(0.0,1.0);
                       iip = isinplace(f,3),
                       noise::NoiseProcess{NoiseClass,inplace,F3}= iip ? INPLACE_WHITE_NOISE : WHITE_NOISE)
  SDETestProblem{typeof(u0),eltype(tspan),iip,inplace,NoiseClass,typeof(f),typeof(g),F3}(f,g,u0,analytic,tspan,noise)
end
