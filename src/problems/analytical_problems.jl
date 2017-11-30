struct AnalyticalProblem{uType,tType,isinplace,F,C} <: AbstractAnalyticalProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  callback::C
  function AnalyticalProblem{iip}(f,u0,tspan;
           callback = nothing) where {iip}
    new{typeof(u0),promote_type(map(typeof,tspan)...),iip,
        typeof(f),typeof(callback)}(f,u0,tspan,callback)
  end
end

function AnalyticalProblem(f,u0,tspan;kwargs...)
  iip = DiffEqBase.isinplace(f,6)
  AnalyticalProblem{iip}(f,u0,tspan;kwargs...)
end

export AnalyticalProblem, AbstractAnalyticalProblem
