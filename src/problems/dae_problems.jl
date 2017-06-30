# f(t,u,du,res) = 0
type DAEProblem{uType,duType,tType,isinplace,F,C,D} <: AbstractDAEProblem{uType,duType,tType,isinplace}
  f::F
  u0::uType
  du0::duType
  tspan::Tuple{tType,tType}
  callback::C
  differential_vars::D
end

function DAEProblem(f,u0,du0,tspan; iip = isinplace(f,4),
                    callback = nothing,
                    differential_vars = nothing)
  DAEProblem{typeof(u0),typeof(du0),promote_type(map(typeof,tspan)...),
             iip,typeof(f),typeof(callback),typeof(differential_vars)}(
             f,u0,du0,tspan,callback,differential_vars)
end
