type ConstantLagDDEProblem{uType,tType,lType,isinplace,F,H,C,MM} <: AbstractConstantLagDDEProblem{uType,tType,lType,isinplace}
  f::F
  h::H
  u0::uType
  lags::lType
  tspan::Tuple{tType,tType}
  callback::C
  mass_matrix::MM
end

function ConstantLagDDEProblem(f,h,u0,lags,tspan; iip = isinplace(f,4), callback = nothing,mass_matrix=I)
  ConstantLagDDEProblem{typeof(u0),eltype(tspan),typeof(lags),iip,typeof(f),typeof(h),typeof(callback),typeof(mass_matrix)}(f,h,u0,lags,tspan,callback,mass_matrix)
end

type DDEProblem{uType,tType,lType,isinplace,F,H,C,MM} <: AbstractDDEProblem{uType,tType,lType,isinplace}
  f::F
  h::H
  u0::uType
  lags::lType
  tspan::Tuple{tType,tType}
  callback::C
  mass_matrix::MM
end

function DDEProblem(f,h,u0,lags,tspan; iip = isinplace(f,4), callback = nothing,mass_matrix=I)
  DDEProblem{typeof(u0),eltype(tspan),typeof(lags),iip,typeof(f),typeof(h),typeof(callback),typeof(mass_matrix)}(f,h,u0,lags,tspan,callback,mass_matrix)
end
