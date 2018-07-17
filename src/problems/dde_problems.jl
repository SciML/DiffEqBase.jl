struct DDEProblem{uType,tType,lType,lType2,isinplace,P,F,H,C,MM} <:
                          AbstractDDEProblem{uType,tType,lType,isinplace}
  f::F
  u0::uType
  h::H
  tspan::tType
  p::P
  constant_lags::lType
  dependent_lags::lType2
  callback::C
  mass_matrix::MM
  neutral::Bool
  @add_kwonly function DDEProblem(f::AbstractDDEFunction,
                                 u0,h,tspan,p=nothing;
                                 constant_lags=[],
                                 dependent_lags=[],
                                 mass_matrix=I,
                                 neutral = mass_matrix == I ?
                                           false : det(mass_matrix)!=1,
                                 callback = nothing)
    _tspan = promote_tspan(tspan)
    new{typeof(u0),typeof(_tspan),
               typeof(constant_lags),typeof(dependent_lags),
               isinplace(f),typeof(p),
               typeof(f),typeof(h),typeof(callback),
               typeof(mass_matrix)}(f,u0,h,_tspan,p,
                                    constant_lags,
                                    dependent_lags,callback,
                                    mass_matrix,neutral)
  end

  function DDEProblem{iip}(f,u0,h,tspan,p=nothing;kwargs...) where {iip}
    DDEProblem(convert(DDEFunction{iip},f),u0,h,tspan,p;kwargs...)
  end

end

function DDEProblem(f,u0,h,tspan,p=nothing;kwargs...)
  DDEProblem(convert(DDEFunction,f),u0,h,tspan,p;kwargs...)
end
