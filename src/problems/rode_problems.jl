type RODEProblem{uType,tType,isinplace,F,C} <: AbstractRODEProblem{uType,tType,isinplace,F}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
end

function RODEProblem(f,u0,tspan; iip = isinplace(f,3),callback=CallbackSet())
  RODEProblem{typeof(u0),eltype(tspan),iip,typeof(f),typeof(calback)}(f,u0,tspan,callback)
end
