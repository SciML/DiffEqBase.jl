struct AnalyticalProblem{uType,tType,isinplace,P,F,C} <: AbstractAnalyticalProblem{uType,tType,isinplace}
  f::F
  u0::uType
  tspan::tType
  p::P
  callback::C
  @add_kwonly function AnalyticalProblem{iip}(f,u0,tspan,p=NullParameters();
           callback = nothing) where {iip}
    _tspan = promote_tspan(tspan)
    new{typeof(u0),typeof(_tspan),iip,typeof(p),
        typeof(f),typeof(callback)}(f,u0,_tspan,p,callback)
  end
end

function AnalyticalProblem(f,u0,tspan,p=NullParameters();kwargs...)
  iip = DiffEqBase.isinplace(f,4)
  AnalyticalProblem{iip}(f,u0,tspan,p;kwargs...)
end

export AnalyticalProblem, AbstractAnalyticalProblem
