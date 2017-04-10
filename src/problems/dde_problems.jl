type ConstantLagDDEProblem{uType,tType,lType,isinplace,F,H,C} <: AbstractConstantLagDDEProblem{uType,tType,lType,isinplace,F,H}
  f::F
  h::H
  u0::uType
  lags::lType
  tspan::Tuple{tType,tType}
  callback::C
end

function ConstantLagDDEProblem(f,h,u0,lags,tspan; iip = isinplace(f,4), callback = CallbackSet())
  ConstantLagDDEProblem{typeof(u0),eltype(tspan),typeof(lags),iip,typeof(f),typeof(h),typeof(callback)}(f,h,u0,lags,tspan,callback)
end

type DDEProblem{uType,tType,lType,isinplace,F,H,C} <: AbstractDDEProblem{uType,tType,lType,isinplace,F,H}
  f::F
  h::H
  u0::uType
  lags::lType
  tspan::Tuple{tType,tType}
  callback::C
end

function DDEProblem(f,h,u0,lags,tspan; iip = isinplace(f,4), callback = CallbackSet())
  DDEProblem{typeof(u0),eltype(tspan),typeof(lags),iip,typeof(f),typeof(h),typeof(callback)}(f,h,u0,lags,tspan,callback)
end
