type ODESolution{uType,uType2,DType,tType,rateType,P,A,IType} <: AbstractODESolution
  u::uType
  u_analytic::uType2
  errors::DType
  t::tType
  k::rateType
  prob::P
  alg::A
  interp::IType
  dense::Bool
  tslocation::Int
  retcode::Symbol
end
(sol::ODESolution)(t,deriv::Type=Val{0};idxs=nothing) = sol.interp(t,idxs,deriv)
(sol::ODESolution)(v,t,deriv::Type=Val{0};idxs=nothing) = sol.interp(v,t,idxs,deriv)

function build_solution{uType,tType,isinplace}(
        prob::Union{AbstractODEProblem{uType,tType,isinplace},AbstractDDEProblem{uType,tType,isinplace}},
        alg,t,u;dense=false,timeseries_errors=true,dense_errors=true,
        calculate_error = true,
        k=[],interp = (tvals) -> nothing, retcode = :Default, kwargs...)

  if has_analytic(prob.f)
    u_analytic = Vector{uType}(0)
    errors = Dict{Symbol,eltype(prob.u0)}()
    sol = ODESolution(u,u_analytic,errors,t,k,prob,alg,interp,dense,0,retcode)
    if calculate_error
      calculate_solution_errors!(sol;timeseries_errors=timeseries_errors,dense_errors=dense_errors)
    end
    return sol
  else
    return ODESolution(u,nothing,nothing,t,k,prob,alg,interp,dense,0,retcode)
  end
end

function calculate_solution_errors!(sol::AbstractODESolution;fill_uanalytic=true,timeseries_errors=true,dense_errors=true)
  if fill_uanalytic
    for i in 1:size(sol.u,1)
      push!(sol.u_analytic,sol.prob.f(Val{:analytic},sol.t[i],sol.prob.u0))
    end
  end

  save_everystep = length(sol.u) > 2
  if !isempty(sol.u_analytic)
    sol.errors[:final] = mean(abs.(sol.u[end]-sol.u_analytic[end]))

    if save_everystep && timeseries_errors
      sol.errors[:l∞] = maximum(vecvecapply((x)->abs.(x),sol.u-sol.u_analytic))
      sol.errors[:l2] = sqrt(mean(vecvecapply((x)->float.(x).^2,sol.u-sol.u_analytic)))
      if sol.dense && dense_errors
        densetimes = collect(linspace(sol.t[1],sol.t[end],100))
        interp_u = sol(densetimes)
        interp_analytic = [sol.prob.f(Val{:analytic},t,sol.u[1]) for t in densetimes]
        sol.errors[:L∞] = maximum(vecvecapply((x)->abs.(x),interp_u-interp_analytic))
        sol.errors[:L2] = sqrt(mean(vecvecapply((x)->float.(x).^2,interp_u-interp_analytic)))
      end
    end
  end
end

function build_solution(sol::AbstractODESolution,u_analytic,errors)
  ODESolution(sol.u,u_analytic,errors,sol.t,sol.k,sol.prob,sol.alg,sol.interp,sol.dense,sol.tslocation,sol.retcode)
end
