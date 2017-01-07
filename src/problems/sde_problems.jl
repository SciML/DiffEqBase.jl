type SDEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3} <: AbstractSDEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3}
  f::F
  g::F2
  u0::uType
  tspan::Tuple{tType,tType}
  noise::NoiseProcess{NoiseClass,F3}
end

type SDETestProblem{uType,tType,isinplace,NoiseClass,F,F2,F3} <: AbstractSDETestProblem{uType,tType,isinplace,NoiseClass,F,F2,F3}
  f::F
  g::F2
  u0::uType
  analytic::Base.Callable
  tspan::Tuple{tType,tType}
  noise::NoiseProcess{NoiseClass,F3}
end

function SDEProblem{NoiseClass,F3}(f,g,u0,tspan,noise::NoiseProcess{NoiseClass,F3}=WHITE_NOISE;iip = isinplace(f,3))
  SDEProblem{typeof(u0),eltype(tspan),iip,NoiseClass,typeof(f),typeof(g),F3}(f,g,u0,tspan,noise)
end

function SDETestProblem{NoiseClass,F3}(f,g,u0,analytic,tspan=(0.0,1.0),noise::NoiseProcess{NoiseClass,F3}=WHITE_NOISE;iip = isinplace(f,3))
  SDETestProblem{typeof(u0),eltype(tspan),iip,NoiseClass,typeof(f),typeof(g),F3}(f,g,u0,analytic,tspan,noise)
end

function print{uType,tType,isinplace,NoiseClass,F,F2,F3}(io::IO, prob::AbstractSDEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3})
  println(io,"AbstractSDEProblem")
  println(io,"Independent Variable Type: $uType")
  println(io,"Depdendent Variable Type: $tType")
  println(io,"Function is in-place? $isinplace")
  nothing
end

function show{uType,tType,isinplace,NoiseClass,F,F2,F3}(io::IO,prob::AbstractSDEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3})
  println(io,"AbstractSDEProblem{$uType,$tType,$isinplace}")
  nothing
end
