struct AnalyticalProblem{uType,auxType,tType,isinplace,F,C} <: AbstractAnalyticalProblem{uType,tType,isinplace}
  f::F
  u0::uType
  aux::auxType
  tspan::Tuple{tType,tType}
  callback::C
  function AnalyticalProblem{iip}(f,u0,aux,tspan;
           callback = nothing) where {iip}
    new{typeof(u0),typeof(aux),promote_type(map(typeof,tspan)...),iip,
        typeof(f),typeof(callback)}(f,u0,aux,tspan,callback)
  end
end

function AnalyticalProblem(f,u0,aux,tspan;kwargs...)
  iip = DiffEqBase.isinplace(f,6)
  AnalyticalProblem{iip}(f,u0,aux,tspan;kwargs...)
end

export AnalyticalProblem, AbstractAnalyticalProblem
