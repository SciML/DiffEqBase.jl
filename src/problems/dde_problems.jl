struct DDEProblem{uType,tType,lType,lType2,isinplace,F,H,C,MM} <:
                          AbstractDDEProblem{uType,tType,lType,isinplace}
  f::F
  h::H
  u0::uType
  constant_lags::lType
  dependent_lags::lType2
  tspan::Tuple{tType,tType}
  callback::C
  mass_matrix::MM
  neutral::Bool
  function DDEProblem{isinplace}(f,h,u0,tspan,constant_lags=nothing,dependent_lags=nothing;
                                 mass_matrix=I,
                                 neutral = mass_matrix == I ?
                                           false : det(mass_matrix)!=1,
                                 callback = nothing) where {isinplace}
    new{typeof(u0),promote_type(map(typeof,tspan)...),
               typeof(constant_lags),typeof(dependent_lags),
               isinplace,typeof(f),typeof(h),typeof(callback),
               typeof(mass_matrix)}(f,h,u0,constant_lags,
                                    dependent_lags,tspan,callback,mass_matrix,neutral)
  end
end

function DDEProblem(f,h,u0,tspan,constant_lags=nothing,dependent_lags=nothing;kwargs...)
  iip = typeof(f)<: Tuple ? isinplace(f[2],4) : isinplace(f,4)
  DDEProblem{iip}(f,h,u0,tspan,constant_lags,dependent_lags;kwargs...)
end
