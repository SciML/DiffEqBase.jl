type DAEProblem{uType,duType,tType,isinplace,F} <: AbstractDAEProblem{uType,duType,tType,isinplace,F}
  f::F
  u0::uType
  du0::duType
  tspan::Tuple{tType,tType}
end

type DAETestProblem{uType,duType,tType,isinplace,F} <: AbstractDAETestProblem{uType,duType,tType,isinplace,F}
  f::F
  u0::uType
  du0::duType
  analytic::Base.Callable
  tspan::Tuple{tType,tType}
end

function DAEProblem(f,u0,du0,tspan)
  iip = isinplace(f,4)
  DAEProblem{typeof(u0),typeof(du0),eltype(tspan),iip,typeof(f)}(f,u0,du0,tspan)
end

function DAETestProblem(f,u0,du0,analytic,tspan=(0.0,1.0))
  iip = isinplace(f,4)
  DAETestProblem{typeof(u0),typeof(du0),eltype(tspan),iip,typeof(f)}(f,u0,du0,analytic,tspan)
end

function print{uType,duType,tType,isinplace,F}(io::IO, prob::AbstractDAEProblem{uType,duType,tType,isinplace,F})
  println(io,"AbstractDAEProblem")
  println(io,"Independent Variable Types: $uType")
  println(io,"Depdendent Variable Type: $tType")
  println(io,"Function is in-place? $isinplace")
  nothing
end

#=
function show{uType,tType,isinplace,F}(io::IO,prob::AbstractODEProblem{uType,tType,isinplace,F})
  println(io,"AbstractDAEProblem{$uType,$tType,$isinplace}")
  nothing
end
=#
