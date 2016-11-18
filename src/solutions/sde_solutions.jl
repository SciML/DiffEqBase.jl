### Concrete Types

type SDESolution{uType,tType,randType,P,A} <: AbstractSDESolution
  u::uType
  t::tType
  W::randType
  prob::P
  alg::A
  maxstacksize::Int
  dense::Bool
  tslocation::Int
end

type SDETestSolution{uType,uType2,uEltype,tType,randType,P,A} <: AbstractSDESolution
  u::uType
  u_analytic::uType2
  errors::Dict{Symbol,uEltype}
  t::tType
  W::randType
  prob::P
  alg::A
  maxstacksize::Int
  dense::Bool
  tslocation::Int
end

function build_solution{uType,tType,isinplace,NoiseClass,F,F2,F3}(
        prob::AbstractSDEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3},
        alg,t,u;W=[],maxstacksize=0,kwargs...)
  SDESolution(u,t,W,prob,alg,maxstacksize,false,0)
end

function build_solution{uType,tType,isinplace,NoiseClass,F,F2,F3}(
        prob::AbstractSDETestProblem{uType,tType,isinplace,NoiseClass,F,F2,F3},
        alg,t,u;W=[],timeseries_errors=true,maxstacksize=0,kwargs...)

  u_analytic = Vector{uType}(0)
  for i in 1:size(u,1)
    push!(u_analytic,prob.analytic(t[i],prob.u0,W[i]))
  end

  save_timeseries = length(u) > 2

  errors = Dict{Symbol,eltype(u[1])}()
  if !isempty(u_analytic)
    errors[:final] = mean(abs.(u[end]-u_analytic[end]))
    if save_timeseries && timeseries_errors
      errors[:l∞] = maximum(vecvecapply((x)->abs.(x),u-u_analytic))
      errors[:l2] = sqrt(mean(vecvecapply((x)->float.(x).^2,u-u_analytic)))
    end
  end
  SDETestSolution(u,u_analytic,errors,t,W,prob,alg,maxstacksize,false,0)
end

function build_solution(sol::AbstractSDESolution,u_analytic,errors)
  SDETestSolution(sol.u,u_analytic,errors,sol.t,sol.W,sol.prob,sol.alg,sol.maxstacksize,sol.dense,sol.tslocation)
end
