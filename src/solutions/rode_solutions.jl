### Concrete Types

type RODESolution{uType,tType,IType,randType,P,A} <: AbstractRODESolution
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
  retcode::Symbol
end
(sol::RODESolution)(t,deriv::Type=Val{0};idxs=nothing) = sol.interp(t,idxs,deriv)
(sol::RODESolution)(v,t,deriv::Type=Val{0};idxs=nothing) = sol.interp(v,t,idxs,deriv)

type RODETestSolution{uType,uType2,uEltype,tType,IType,randType,P,A} <: AbstractRODETestSolution
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
  retcode::Symbol
end
(sol::RODETestSolution)(t,deriv::Type=Val{0};idxs=nothing) = sol.interp(t,idxs,deriv)
(sol::RODETestSolution)(v,t,deriv::Type=Val{0};idxs=nothing) = sol.interp(v,t,idxs,deriv)

function build_solution{uType,tType,isinplace,NoiseClass,F,F2,F3}(
        prob::AbstractRODEProblem{uType,tType,isinplace,NoiseClass,F,F2,F3},
        alg,t,u;W=[],maxstacksize=0,maxstacksize2=0,
        interp = (tvals) -> nothing,
        dense = false,retcode = :Default, kwargs...)
  RODESolution(u,t,W,prob,alg,maxstacksize,maxstacksize2,interp,dense,0,retcode)
end

function build_solution{uType,tType,isinplace,NoiseClass,F,F2,F3}(
        prob::Union{AbstractRODETestProblem{uType,tType,isinplace,NoiseClass,F,F2,F3},AbstractSDETestProblem{uType,tType,isinplace,NoiseClass,F,F2,F3}},
        alg,t,u;W=[],timeseries_errors=true,dense_errors=false,calculate_error=true,
        maxstacksize=0,maxstacksize2=0,
        interp = (tvals) -> nothing,
        dense = false, retcode = :Default, kwargs...)

  u_analytic = Vector{uType}(0)
  save_timeseries = length(u) > 2
  errors = Dict{Symbol,eltype(u[1])}()
  sol = RODETestSolution(u,u_analytic,errors,t,W,prob,alg,maxstacksize,maxstacksize2,interp,dense,0,retcode)
  if calculate_error
    calculate_solution_errors!(sol;timeseries_errors=timeseries_errors,dense_errors=dense_errors)
  end
  sol
end

function calculate_solution_errors!(sol::AbstractRODETestSolution;fill_uanalytic=true,timeseries_errors=true,dense_errors=true)
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

function build_solution(sol::AbstractRODESolution,u_analytic,errors)
  RODETestSolution(sol.u,u_analytic,errors,sol.t,sol.W,sol.prob,sol.alg,sol.maxstacksize,sol.dense,sol.tslocation,sol.retcode)
end
