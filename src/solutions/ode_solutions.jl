### Concrete Types

type ODESolution{uType,tType,rateType,P,A,IType} <: AbstractODESolution
  u::uType
  t::tType
  k::rateType
  prob::P
  alg::A
  interp::IType
  dense::Bool
  tslocation::Int
end
(sol::ODESolution)(t) = sol.interp(t)

type ODETestSolution{uType,uType2,uEltype,tType,rateType,P,A,IType} <: AbstractODETestSolution
  u::uType
  u_analytic::uType2
  errors::Dict{Symbol,uEltype}
  t::tType
  k::rateType
  prob::P
  alg::A
  interp::IType
  dense::Bool
  tslocation::Int
end
(sol::ODETestSolution)(t) = sol.interp(t)

function build_solution{uType,tType,isinplace}(
        prob::AbstractODEProblem{uType,tType,isinplace},
        alg,t,u;dense=false,
        k=[],interp = (tvals) -> nothing,kwargs...)
  ODESolution(u,t,k,prob,alg,interp,dense,0)
end

function build_solution{uType,tType,isinplace}(
        prob::AbstractODETestProblem{uType,tType,isinplace},
        alg,t,u;dense=false,
        k=[],interp = (tvals) -> nothing,
        timeseries_errors=true,dense_errors=true,
        calculate_error = true,kwargs...)
  u_analytic = Vector{uType}(0)
  errors = Dict{Symbol,eltype(u[1])}()
  sol = ODETestSolution(u,u_analytic,errors,t,k,prob,alg,interp,dense,0)
  if calculate_error
    calculate_solution_errors!(sol;timeseries_errors=timeseries_errors,dense_errors=dense_errors)
  end
  sol
end

function calculate_solution_errors!(sol::AbstractODETestSolution;fill_uanalytic=true,timeseries_errors=true,dense_errors=true)
  if fill_uanalytic
    for i in 1:size(sol.u,1)
      push!(sol.u_analytic,sol.prob.analytic(sol.t[i],sol.prob.u0))
    end
  end

  save_timeseries = length(sol.u) > 2
  if !isempty(sol.u_analytic)
    sol.errors[:final] = mean(abs.(sol.u[end]-sol.u_analytic[end]))

    if save_timeseries && timeseries_errors
      sol.errors[:l∞] = maximum(vecvecapply((x)->abs.(x),sol.u-sol.u_analytic))
      sol.errors[:l2] = sqrt(mean(vecvecapply((x)->float.(x).^2,sol.u-sol.u_analytic)))
      if sol.dense && dense_errors
        densetimes = collect(linspace(sol.t[1],sol.t[end],100))
        interp_u = sol(densetimes)
        interp_analytic = [sol.prob.analytic(t,sol.u[1]) for t in densetimes]
        sol.errors[:L∞] = maximum(vecvecapply((x)->abs.(x),interp_u-interp_analytic))
        sol.errors[:L2] = sqrt(mean(vecvecapply((x)->float.(x).^2,interp_u-interp_analytic)))
      end
    end
  end
end

function build_solution(sol::AbstractODESolution,u_analytic,errors)
  ODETestSolution(sol.u,u_analytic,errors,sol.t,sol.k,sol.prob,sol.alg,sol.interp,sol.dense,sol.tslocation)
end
