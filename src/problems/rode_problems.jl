type RODEProblem{uType,tType,isinplace,F} <: AbstractRODEProblem{uType,tType,isinplace,F}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
end

type RODETestProblem{uType,AType,tType,isinplace,F} <: AbstractRODETestProblem{uType,tType,isinplace,F}
  f::F
  u0::uType
  analytic::AType
  tspan::Tuple{tType,tType}
end

function RODEProblem(f,u0,tspan; iip = isinplace(f,3))
  RODEProblem{typeof(u0),eltype(tspan),iip,typeof(f)}(f,u0,tspan)
end

function RODETestProblem(f,u0,analytic,tspan=(0.0,1.0); iip = isinplace(f,3))
  RODETestProblem{typeof(u0),typeof(analytic),eltype(tspan),iip,typeof(f)}(f,u0,analytic,tspan)
end
