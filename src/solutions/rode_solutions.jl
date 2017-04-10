### Concrete Types

type RODESolution{uType,uType2,DType,tType,IType,randType,P,A} <: AbstractRODESolution
  u::uType
  u_analytic::uType2
  errors::DType
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

function build_solution{uType,tType,isinplace,NoiseClass}(
        prob::AbstractRODEProblem{uType,tType,isinplace,NoiseClass},
        alg,t,u;W=[],timeseries_errors=true,dense_errors=false,calculate_error=true,
        maxstacksize=0,maxstacksize2=0,
        interp = (tvals) -> nothing,
        dense = false, retcode = :Default, kwargs...)

  if has_analytic(prob.f)
    u_analytic = Vector{uType}(0)
    errors = Dict{Symbol,eltype(u[1])}()
    sol = RODESolution(u,u_analytic,errors,t,W,prob,alg,maxstacksize,maxstacksize2,interp,dense,0,retcode)
    if calculate_error
      calculate_solution_errors!(sol;timeseries_errors=timeseries_errors,dense_errors=dense_errors)
    end
    return sol
  else
    return RODESolution(u,nothing,nothing,t,W,prob,alg,maxstacksize,maxstacksize2,interp,dense,0,retcode)
  end
end

function calculate_solution_errors!(sol::AbstractRODESolution;fill_uanalytic=true,timeseries_errors=true,dense_errors=true)
  if fill_uanalytic
    for i in 1:size(sol.u,1)
      push!(sol.u_analytic,sol.prob.f(Val{:analytic},sol.t[i],sol.prob.u0,sol.W[i]))
    end
  end

  save_everystep = length(sol.u) > 2

  if !isempty(sol.u_analytic)
    sol.errors[:final] = mean(abs.(sol.u[end]-sol.u_analytic[end]))
    if save_everystep && timeseries_errors
      sol.errors[:l∞] = maximum(vecvecapply((x)->abs.(x),sol.u-sol.u_analytic))
      sol.errors[:l2] = sqrt(mean(vecvecapply((x)->float.(x).^2,sol.u-sol.u_analytic)))
    end
  end
end

function build_solution(sol::AbstractRODESolution,u_analytic,errors)
  RODESolution(sol.u,u_analytic,errors,sol.t,sol.W,sol.prob,sol.alg,sol.maxstacksize,sol.dense,sol.tslocation,sol.retcode)
end
