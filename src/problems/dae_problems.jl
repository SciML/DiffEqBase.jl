# f(t,u,du,res) = 0
struct DAEProblem{uType,duType,tType,isinplace,P,F,C,D} <: AbstractDAEProblem{uType,duType,tType,isinplace}
  f::F
  du0::duType
  u0::uType  
  tspan::Tuple{tType,tType}
  p::P
  callback::C
  differential_vars::D
  function DAEProblem{iip}(f,du0,u0,tspan,p=nothing;
                      callback = nothing,
                      differential_vars = nothing) where {iip}
    new{typeof(u0),typeof(du0),promote_type(map(typeof,tspan)...),
               iip,typeof(p),typeof(f),typeof(callback),
               typeof(differential_vars)}(
               f,du0,u0,tspan,p,callback,differential_vars)
  end
end

function DAEProblem(f,du0,u0,tspan,p=nothing;kwargs...)
  iip = typeof(f)<: Tuple ? isinplace(f[2],5) : isinplace(f,5)
  DAEProblem{iip}(f,du0,u0,tspan,p;kwargs...)
end
