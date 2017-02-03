type ConstantLagDDEProblem{uType,tType,lType,isinplace,F,H} <: AbstractConstantLagDDEProblem{uType,tType,lType,isinplace,F,H}
  f::F
  h::H
  u0::uType
  lags::lType
  tspan::Tuple{tType,tType}
end

type ConstantLagDDETestProblem{uType,tType,lType,isinplace,F,H} <: AbstractConstantLagDDETestProblem{uType,tType,lType,isinplace,F,H}
  f::F
  h::H
  u0::uType
  analytic::Base.Callable
  lags::lType
  tspan::Tuple{tType,tType}
end

function ConstantLagDDEProblem(f,h,u0,lags,tspan; iip = isinplace(f,4))
  ConstantLagDDEProblem{typeof(u0),eltype(tspan),typeof(lags),iip,typeof(f),typeof(h)}(f,h,u0,lags,tspan)
end

function ConstantLagDDETestProblem(f,h,u0,analytic,lags,tspan=(0.0,1.0); iip = isinplace(f,4))
  ConstantLagDDETestProblem{typeof(u0),eltype(tspan),typeof(lags),iip,typeof(f),typeof(h)}(f,h,u0,analytic,lags,tspan)
end

type DDEProblem{uType,tType,lType,isinplace,F,H} <: AbstractDDEProblem{uType,tType,lType,isinplace,F,H}
  f::F
  h::H
  u0::uType
  lags::lType
  tspan::Tuple{tType,tType}
end

type DDETestProblem{uType,tType,lType,isinplace,F,H} <: AbstractDDETestProblem{uType,tType,lType,isinplace,F,H}
  f::F
  h::H
  u0::uType
  analytic::Base.Callable
  lags::lType
  tspan::Tuple{tType,tType}
end

function DDEProblem(f,h,u0,lags,tspan; iip = isinplace(f,4))
  DDEProblem{typeof(u0),eltype(tspan),typeof(lags),iip,typeof(f),typeof(h)}(f,h,u0,lags,tspan)
end

function DDETestProblem(f,h,u0,analytic,lags,tspan=(0.0,1.0); iip = isinplace(f,4))
  DDETestProblem{typeof(u0),eltype(tspan),typeof(lags),iip,typeof(f),typeof(h)}(f,h,u0,analytic,lags,tspan)
end

#=
function print{uType,tType,isinplace,F,H}(io::IO, prob::AbstractDDEProblem{uType,tType,isinplace,F,H})
  println(io,"AbstractDDEProblem")
  println(io,"Independent Variable Type: $uType")
  println(io,"Depdendent Variable Type: $tType")
  println(io,"Function is in-place? $isinplace")
  nothing
end
=#

#=
function show{uType,tType,isinplace,F}(io::IO,prob::AbstractODEProblem{uType,tType,isinplace,F})
  println(io,"AbstractODEProblem{$uType,$tType,$isinplace}")
  nothing
end
=#
