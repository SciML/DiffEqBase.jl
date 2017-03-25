type IMEXProblem{uType,tType,isinplace,F,G} <: AbstractODEProblem{uType,tType,isinplace,F,G}
  f::F
  g::G
  u0::uType
  tspan::Tuple{tType,tType}
end

type IMEXTestProblem{uType,AType,tType,isinplace,F,G} <: AbstractODETestProblem{uType,tType,isinplace,F,G}
  f::F
  g::G
  u0::uType
  analytic::AType
  tspan::Tuple{tType,tType}
end

function IMEXProblem(f,g,u0,tspan; iip = isinplace(f,3))
  IMEXProblem{typeof(u0),eltype(tspan),iip,typeof(f),typeof(g)}(f,g,u0,tspan)
end

function IMEXTestProblem(f,g,u0,analytic,tspan=(0.0,1.0); iip = isinplace(f,3))
  IMEXTestProblem{typeof(u0),typeof(analytic),eltype(tspan),iip,typeof(f),typeof(g)}(f,g,u0,analytic,tspan)
end

#=
function print{uType,tType,isinplace,F,G}(io::IO, prob::AbstractODEProblem{uType,tType,isinplace,F,G})
  println(io,"AbstractIMEXProblem")
  println(io,"Independent Variable Type: $uType")
  println(io,"Depdendent Variable Type: $tType")
  println(io,"Function is in-place? $isinplace")
  nothing
end
=#

#=
function show{uType,tType,isinplace,F,G}(io::IO,prob::AbstractODEProblem{uType,tType,isinplace,F,G})
  println(io,"AbstractIMEXProblem{$uType,$tType,$isinplace}")
  nothing
end
=#