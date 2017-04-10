type ODEProblem{uType,tType,isinplace,F,C} <: AbstractODEProblem{uType,tType,isinplace,F}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
end

function ODEProblem(f,u0,tspan; iip = isinplace(f,3),callback=CallbackSet())
  ODEProblem{typeof(u0),eltype(tspan),iip,typeof(f),typeof(callback)}(f,u0,tspan,callback)
end
