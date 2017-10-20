# f(t,u,du,res) = 0
struct DAEProblem{uType,duType,tType,isinplace,F,C,D} <: AbstractDAEProblem{uType,duType,tType,isinplace}
  f::F
  u0::uType
  du0::duType
  tspan::Tuple{tType,tType}
  callback::C
  differential_vars::D
  function DAEProblem{iip}(f,u0,du0,tspan;
                      callback = nothing,
                      differential_vars = nothing) where {iip}
    new{typeof(u0),typeof(du0),promote_type(map(typeof,tspan)...),
               iip,typeof(f),typeof(callback),
               typeof(differential_vars)}(
               f,u0,du0,tspan,callback,differential_vars)
  end
end

function DAEProblem(f,u0,du0,tspan;kwargs...)
  iip = typeof(f)<: Tuple ? isinplace(f[2],4) : isinplace(f,4)
  DAEProblem{iip}(f,u0,du0,tspan;kwargs...)
end
