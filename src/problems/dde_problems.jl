type ConstantLagDDEProblem{uType,tType,lType,isinplace,F,H,C,MM} <: AbstractConstantLagDDEProblem{uType,tType,lType,isinplace}
  f::F
  h::H
  u0::uType
  lags::lType
  tspan::Tuple{tType,tType}
  callback::C
  mm::MM
end

function ConstantLagDDEProblem(f,h,u0,lags,tspan; iip = isinplace(f,4), callback = nothing,mm=I)
  ConstantLagDDEProblem{typeof(u0),eltype(tspan),typeof(lags),iip,typeof(f),typeof(h),typeof(callback),typeof(mm)}(f,h,u0,lags,tspan,callback,mm)
end

type DDEProblem{uType,tType,lType,isinplace,F,H,C,MM} <: AbstractDDEProblem{uType,tType,lType,isinplace}
  f::F
  h::H
  u0::uType
  lags::lType
  tspan::Tuple{tType,tType}
  callback::C
  mm::MM
end

function DDEProblem(f,h,u0,lags,tspan; iip = isinplace(f,4), callback = nothing,mm=I)
  DDEProblem{typeof(u0),eltype(tspan),typeof(lags),iip,typeof(f),typeof(h),typeof(callback),typeof(mm)}(f,h,u0,lags,tspan,callback,mm)
end
