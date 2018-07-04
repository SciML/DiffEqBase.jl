struct DDEProblem{uType,tType,lType,lType2,isinplace,P,J,F,H,C,MM} <:
                          AbstractDDEProblem{uType,tType,lType,isinplace}
  f::F
  u0::uType
  h::H
  tspan::Tuple{tType,tType}
  p::P
  jac_prototype::J
  constant_lags::lType
  dependent_lags::lType2
  callback::C
  mass_matrix::MM
  neutral::Bool
  @add_kwonly function DDEProblem{isinplace}(f,u0,h,tspan,p=nothing;
                                 jac_prototype = nothing,
                                 constant_lags=[],
                                 dependent_lags=[],
                                 mass_matrix=I,
                                 neutral = mass_matrix == I ?
                                           false : det(mass_matrix)!=1,
                                 callback = nothing) where {isinplace}
    _tspan = promote_tspan(tspan)
    new{typeof(u0),typeof(_tspan),
               typeof(constant_lags),typeof(dependent_lags),
               isinplace,typeof(p),typeof(jac_prototype),
               typeof(f),typeof(h),typeof(callback),
               typeof(mass_matrix)}(f,u0,h,_tspan,p,
                                    jac_prototype,
                                    constant_lags,
                                    dependent_lags,callback,
                                    mass_matrix,neutral)
  end
end

function DDEProblem(f,u0,h,tspan,p=nothing;kwargs...)
  iip = typeof(f)<: Tuple ? isinplace(f[2],5) : isinplace(f,5)
  DDEProblem{iip}(f,u0,h,tspan,p;kwargs...)
end
