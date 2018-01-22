struct DDEProblem{uType,tType,lType,lType2,isinplace,P,F,H,C,MM} <:
                          AbstractDDEProblem{uType,tType,lType,isinplace}
  f::F
  h::H
  u0::uType
  constant_lags::lType
  dependent_lags::lType2
  tspan::Tuple{tType,tType}
  p::P
  callback::C
  mass_matrix::MM
  neutral::Bool
  function DDEProblem{isinplace}(f,h,u0,tspan,p=nothing,
                                 constant_lags=nothing,
                                 dependent_lags=nothing;
                                 mass_matrix=I,
                                 neutral = mass_matrix == I ?
                                           false : det(mass_matrix)!=1,
                                 callback = nothing) where {isinplace}
    new{typeof(u0),promote_type(map(typeof,tspan)...),
               typeof(constant_lags),typeof(dependent_lags),
               isinplace,typeof(p),typeof(f),typeof(h),typeof(callback),
               typeof(mass_matrix)}(f,h,u0,constant_lags,
                                    dependent_lags,tspan,p,callback,
                                    mass_matrix,neutral)
  end
end

function DDEProblem(f,h,u0,tspan,p=nothing,
                    constant_lags=nothing,
                    dependent_lags=nothing;kwargs...)
  iip = typeof(f)<: Tuple ? isinplace(f[2],5) : isinplace(f,5)
  DDEProblem{iip}(f,h,u0,tspan,p,constant_lags,dependent_lags;kwargs...)
end
