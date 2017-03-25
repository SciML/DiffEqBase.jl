type SplitODEProblem{uType,tType,isinplace,F} <: AbstractSplitODEProblem{uType,tType,isinplace,F}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
end

type SplitODETestProblem{uType,AType,tType,isinplace,F} <: AbstractSplitODETestProblem{uType,tType,isinplace,F}
  f::F
  u0::uType
  analytic::AType
  tspan::Tuple{tType,tType}
end

function SplitODEProblem(f,u0,tspan; iip = isinplace(f[2],3))
  SplitODEProblem{typeof(u0),eltype(tspan),iip,typeof(f)}(f,u0,tspan)
end

function SplitODETestProblem(f,u0,analytic,tspan=(0.0,1.0); iip = isinplace(f[2],3))
  SplitODETestProblem{typeof(u0),typeof(analytic),eltype(tspan),iip,typeof(f)}(f,u0,analytic,tspan)
end

#=
function print{uType,tType,isinplace,F}(io::IO, prob::AbstractSplitODEProblem{uType,tType,isinplace,F})
  println(io,"AbstractSplitODEProblem")
  println(io,"Independent Variable Type: $uType")
  println(io,"Depdendent Variable Type: $tType")
  println(io,"Function is in-place? $isinplace")
  nothing
end
=#

#=
function show{uType,tType,isinplace,F}(io::IO,prob::AbstractSplitODEProblem{uType,tType,isinplace,F})
  println(io,"AbstractSplitODEProblem{$uType,$tType,$isinplace}")
  nothing
end
=#
