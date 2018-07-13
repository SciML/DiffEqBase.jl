# f(t,u,du,res) = 0
struct DAEProblem{uType,duType,tType,isinplace,P,J,F,C,D} <: AbstractDAEProblem{uType,duType,tType,isinplace}
  f::F
  du0::duType
  u0::uType
  tspan::tType
  p::P
  jac_prototype::J
  callback::C
  differential_vars::D
  @add_kwonly function DAEProblem(f::AbstractODEFunction,du0,u0,tspan,p=nothing;
                      jac_prototype = nothing,
                      callback = nothing,
                      differential_vars = nothing)
    _tspan = promote_tspan(tspan)
    new{typeof(u0),typeof(du0),typeof(_tspan),
               iip,typeof(p),typeof(jac_prototype),
               typeof(f),typeof(callback),
               typeof(differential_vars)}(
               f,du0,u0,_tspan,p,jac_prototype,
               callback,differential_vars)
  end
  function DAEProblem{iip}(f,du0,u0,tspan,p=nothing;kwargs...) where {iip}
    DAEProblem(convert(DAEFunction{iip},f),du0,u0,tspan,p;kwargs...)
  end
end

function DAEProblem(f,du0,u0,tspan,p=nothing;kwargs...)
  DAEProblem(convert(DAEFunction,f),du0,u0,tspan,p;kwargs...)
end
