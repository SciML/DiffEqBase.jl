### Concrete Types

type ODESolution{uType,tType,rateType,P,A} <: AbstractODESolution
  u::uType
  t::tType
  k::rateType
  prob::P
  alg::A
  interp::Function
  dense::Bool
  tslocation::Int
end
(sol::ODESolution)(t) = sol.interp(t)

type ODETestSolution{uType,uType2,uEltype,tType,rateType,P,A} <: AbstractODESolution
  u::uType
  u_analytic::uType2
  errors::Dict{Symbol,uEltype}
  t::tType
  k::rateType
  prob::P
  alg::A
  interp::Function
  dense::Bool
  tslocation::Int
end
(sol::ODETestSolution)(t) = sol.interp(t)

function build_ode_solution{uType,tType,isinplace}(
        prob::AbstractODEProblem{uType,tType,isinplace},
        alg,t,u;dense=false,
        k=[],interp = (tvals) -> nothing,kwargs...)
  ODESolution(u,t,k,prob,alg,interp,dense,0)
end

function build_ode_solution{uType,tType,isinplace}(
        prob::AbstractODETestProblem{uType,tType,isinplace},
        alg,t,u;dense=false,
        k=[],interp = (tvals) -> nothing,
        timeseries_errors=true,dense_errors=true,kwargs...)

  u_analytic = Vector{uType}(0)
  for i in 1:size(u,1)
    push!(u_analytic,prob.analytic(t[i],prob.u0))
  end

  save_timeseries = length(u) > 2

  errors = Dict{Symbol,eltype(u[1])}()
  if !isempty(u_analytic)
    errors[:final] = mean(abs.(u[end]-u_analytic[end]))

    if save_timeseries && timeseries_errors
      errors[:l∞] = maximum(vecvecapply((x)->abs.(x),u-u_analytic))
      errors[:l2] = sqrt(mean(vecvecapply((x)->float.(x).^2,u-u_analytic)))
      if dense && dense_errors
        densetimes = collect(linspace(t[1],t[end],100))
        interp_u = interp(densetimes)
        interp_analytic = [prob.analytic(t,u[1]) for t in densetimes]
        errors[:L∞] = maximum(vecvecapply((x)->abs.(x),interp_u-interp_analytic))
        errors[:L2] = sqrt(mean(vecvecapply((x)->float.(x).^2,interp_u-interp_analytic)))
      end
    end
  end
  ODETestSolution(u,u_analytic,errors,t,k,prob,alg,interp,dense,0)
end

function build_ode_solution(sol::AbstractODESolution,u_analytic,errors)
  ODETestSolution(sol.u,u_analytic,errors,sol.t,sol.k,sol.prob,sol.alg,sol.interp,sol.dense,sol.tslocation)
end
