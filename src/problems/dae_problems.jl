type DAEProblem{uType,duType,tType,isinplace,F,C} <: AbstractDAEProblem{uType,duType,tType,isinplace}
  f::F
  u0::uType
  du0::duType
  tspan::Tuple{tType,tType}
  callback::C
end

function DAEProblem(f,u0,du0,tspan; iip = isinplace(f,4), callback = CallbackSet())
  DAEProblem{typeof(u0),typeof(du0),eltype(tspan),iip,typeof(f),typeof(callback)}(f,u0,du0,tspan,callback)
end
