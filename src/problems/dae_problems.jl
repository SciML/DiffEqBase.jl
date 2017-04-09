type DAEProblem{uType,duType,tType,isinplace,F} <: AbstractDAEProblem{uType,duType,tType,isinplace,F}
  f::F
  u0::uType
  du0::duType
  tspan::Tuple{tType,tType}
end

type DAETestProblem{uType,duType,tType,isinplace,F} <: AbstractDAETestProblem{uType,duType,tType,isinplace,F}
  f::F
  u0::uType
  du0::duType
  analytic::Base.Callable
  tspan::Tuple{tType,tType}
end

function DAEProblem(f,u0,du0,tspan; iip = isinplace(f,4))
  DAEProblem{typeof(u0),typeof(du0),eltype(tspan),iip,typeof(f)}(f,u0,du0,tspan)
end

function DAETestProblem(f,u0,du0,analytic,tspan=(0.0,1.0); iip = isinplace(f,4))
  DAETestProblem{typeof(u0),typeof(du0),eltype(tspan),iip,typeof(f)}(f,u0,du0,analytic,tspan)
end
