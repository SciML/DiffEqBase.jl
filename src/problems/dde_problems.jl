type ConstantLagDDEProblem{uType,tType,lType,isinplace,F,H,C,MM} <: AbstractConstantLagDDEProblem{uType,tType,lType,isinplace}
  f::F
  h::H
  u0::uType
  lags::lType
  tspan::Tuple{tType,tType}
  callback::C
  mass_matrix::MM
  function ConstantLagDDEProblem{isinplace}(f,h,u0,lags,tspan;
    callback = nothing,mass_matrix=I) where {isinplace}
    new{typeof(u0),promote_type(map(typeof,tspan)...),
                   typeof(lags),isinplace,typeof(f),typeof(h),
                   typeof(callback),typeof(mass_matrix)}(
                   f,h,u0,lags,tspan,callback,mass_matrix)
  end
end

function ConstantLagDDEProblem(f,h,u0,lags,tspan;kwargs...)
  iip = typeof(f)<: Tuple ? isinplace(f[2],4) : isinplace(f,4)
  ConstantLagDDEProblem{iip}(f,h,u0,lags,tspan;kwargs...)
end

type DDEProblem{uType,tType,lType,isinplace,F,H,C,MM} <: AbstractDDEProblem{uType,tType,lType,isinplace}
  f::F
  h::H
  u0::uType
  lags::lType
  tspan::Tuple{tType,tType}
  callback::C
  mass_matrix::MM
  function DDEProblem{isinplace}(f,h,u0,lags,tspan;
                                 callback = nothing,mass_matrix=I) where {isinplace}
    new{typeof(u0),promote_type(map(typeof,tspan)...),
               typeof(lags),isinplace,typeof(f),typeof(h),typeof(callback),
               typeof(mass_matrix)}(f,h,u0,lags,tspan,callback,mass_matrix)
  end
end

function DDEProblem(f,h,u0,lags,tspan;kwargs...)
  iip = typeof(f)<: Tuple ? isinplace(f[2],4) : isinplace(f,4)
  DDEProblem{iip}(f,h,u0,lags,tspan;kwargs...)
end
