### Concrete Types

type SDESolution{uType,tType,IType,randType,P,A} <: AbstractSDESolution
  u::uType
  t::tType
  W::randType
  prob::P
  alg::A
  maxstacksize::Int
  maxstacksize2::Int
  interp::IType
  dense::Bool
  tslocation::Int
end
(sol::SDESolution)(t,deriv::Type=Val{0};idxs=nothing) = sol.interp(t,idxs,deriv)
(sol::SDESolution)(v,t,deriv::Type=Val{0};idxs=nothing) = sol.interp(v,t,idxs,deriv)

type SDETestSolution{uType,uType2,uEltype,tType,IType,randType,P,A} <: AbstractSDETestSolution
  u::uType
  u_analytic::uType2
  errors::Dict{Symbol,uEltype}
  t::tType
  W::randType
  prob::P
  alg::A
  maxstacksize::Int
  maxstacksize2::Int
  interp::IType
  dense::Bool
  tslocation::Int
end
(sol::SDETestSolution)(t,deriv::Type=Val{0};idxs=nothing) = sol.interp(t,idxs,deriv)
(sol::SDETestSolution)(v,t,deriv::Type=Val{0};idxs=nothing) = sol.interp(v,t,idxs,deriv)

function build_solution{uType,tType,isinplace,NoiseClass,F,F2,F3}(
        prob::AbstractSDEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3},
        alg,t,u;W=[],maxstacksize=0,maxstacksize2=0,
        interp = (tvals) -> nothing,
        dense = false,kwargs...)
  SDESolution(u,t,W,prob,alg,maxstacksize,maxstacksize2,interp,dense,0)
end

function build_solution{uType,tType,isinplace,NoiseClass,F,F2,F3}(
        prob::AbstractSDETestProblem{uType,tType,isinplace,NoiseClass,F,F2,F3},
        alg,t,u;W=[],timeseries_errors=true,dense_errors=false,calculate_error=true,
        maxstacksize=0,maxstacksize2=0,
        interp = (tvals) -> nothing,
        dense = false,kwargs...)

  u_analytic = Vector{uType}(0)
  save_timeseries = length(u) > 2
  errors = Dict{Symbol,eltype(u[1])}()
  sol = SDETestSolution(u,u_analytic,errors,t,W,prob,alg,maxstacksize,maxstacksize2,interp,dense,0)
  if calculate_error
    calculate_solution_errors!(sol;timeseries_errors=timeseries_errors,dense_errors=dense_errors)
  end
  sol
end

function calculate_solution_errors!(sol::AbstractSDETestSolution;fill_uanalytic=true,timeseries_errors=true,dense_errors=true)
  if fill_uanalytic
    for i in 1:size(sol.u,1)
      push!(sol.u_analytic,sol.prob.analytic(sol.t[i],sol.prob.u0,sol.W[i]))
    end
  end

  save_timeseries = length(sol.u) > 2

  if !isempty(sol.u_analytic)
    sol.errors[:final] = mean(abs.(sol.u[end]-sol.u_analytic[end]))
    if save_timeseries && timeseries_errors
      sol.errors[:l∞] = maximum(vecvecapply((x)->abs.(x),sol.u-sol.u_analytic))
      sol.errors[:l2] = sqrt(mean(vecvecapply((x)->float.(x).^2,sol.u-sol.u_analytic)))
    end
  end
end

function build_solution(sol::AbstractSDESolution,u_analytic,errors)
  SDETestSolution(sol.u,u_analytic,errors,sol.t,sol.W,sol.prob,sol.alg,sol.maxstacksize,sol.dense,sol.tslocation)
end
