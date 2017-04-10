type SplitODEProblem{uType,tType,isinplace,F,C} <: AbstractSplitODEProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
end

function SplitODEProblem(f,u0,tspan; iip = isinplace(f[2],3),callback=CallbackSet())
  SplitODEProblem{typeof(u0),eltype(tspan),iip,typeof(f),typeof(callback)}(f,u0,tspan,callback)
end
