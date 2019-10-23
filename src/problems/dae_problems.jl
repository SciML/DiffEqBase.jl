# f(t,u,du,res) = 0
"""
$(TYPEDEF)

TODO
"""
struct DAEProblem{uType,duType,tType,isinplace,P,F,K,D} <: AbstractDAEProblem{uType,duType,tType,isinplace}
  f::F
  du0::duType
  u0::uType
  tspan::tType
  p::P
  kwargs::K
  differential_vars::D
  @add_kwonly function DAEProblem{iip}(f::AbstractDAEFunction{iip},
                      du0,u0,tspan,p=NullParameters();
                      differential_vars = nothing,
                      kwargs...) where {iip}
    _tspan = promote_tspan(tspan)
    new{typeof(u0),typeof(du0),typeof(_tspan),
               isinplace(f),typeof(p),
               typeof(f),typeof(kwargs),
               typeof(differential_vars)}(
               f,du0,u0,_tspan,p,
               kwargs,differential_vars)
  end

  function DAEProblem{iip}(f,du0,u0,tspan,p=NullParameters();kwargs...) where {iip}
    DAEProblem(convert(DAEFunction{iip},f),du0,u0,tspan,p;kwargs...)
  end
end

function DAEProblem(f::AbstractDAEFunction,du0,u0,tspan,p=NullParameters();kwargs...)
  DAEProblem{isinplace(f)}(f,du0,u0,tspan,p;kwargs...)
end

function DAEProblem(f,du0,u0,tspan,p=NullParameters();kwargs...)
  DAEProblem(convert(DAEFunction,f),du0,u0,tspan,p;kwargs...)
end
