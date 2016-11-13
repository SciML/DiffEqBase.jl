type ODEProblem{uType,tType,isinplace,F} <: AbstractODEProblem{uType,tType,isinplace,F}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
end

type ODETestProblem{uType,tType,isinplace,F} <: AbstractODETestProblem{uType,tType,isinplace,F}
  f::F
  u0::uType
  analytic::Base.Callable
  tspan::Tuple{tType,tType}
end

function ODEProblem(f::Base.Callable,u0,tspan)
  isinplace = numparameters(f)>=3
  ODEProblem{typeof(u0),eltype(tspan),isinplace,typeof(f)}(f,u0,tspan)
end

function ODETestProblem(f::Base.Callable,u0,analytic,tspan=(0.0,1.0))
  isinplace = numparameters(f)>=3
  ODETestProblem{typeof(u0),eltype(tspan),isinplace,typeof(f)}(f,u0,analytic,tspan)
end

function print{uType,tType,isinplace,F}(io::IO, prob::AbstractODEProblem{uType,tType,isinplace,F})
  println(io,"AbstractODEProblem")
  println(io,"Independent Variable Type: $uType")
  println(io,"Depdendent Variable Type: $tType")
  println(io,"Function is in-place? $isinplace")
  nothing
end

function show{uType,tType,isinplace,F}(io::IO,prob::AbstractODEProblem{uType,tType,isinplace,F})
  println(io,"AbstractODEProblem{$uType,$tType,$isinplace}")
  nothing
end
