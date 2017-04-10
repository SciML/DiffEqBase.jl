type DAESolution{uType,duType,uType2,DType,tType,ID,P,A} <: AbstractODESolution
  u::uType
  du::duType
  u_analytic::uType2
  errors::DType
  interp_data::ID
  t::tType
  prob::P
  alg::A
  interp::Function
  dense::Bool
  tslocation::Int
  retcode::Symbol
end
(sol::DAESolution)(t,deriv::Type=Val{0};idxs=nothing) = sol.interp(t,idxs,deriv)
(sol::DAESolution)(v,t,deriv::Type=Val{0};idxs=nothing) = sol.interp(v,t,idxs,deriv)

function build_solution{uType,duType,tType,isinplace}(
        prob::AbstractDAEProblem{uType,duType,tType,isinplace},
        alg,t,u;dense=false,du=[],
        interp_data=[],interp = (tvals) -> nothing,
        timeseries_errors=true,dense_errors=true, retcode = :Default, kwargs...)

  if has_analytic(prob.f)
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
    DAESolution(u,du,u_analytic,errors,t,interp_data,prob,alg,interp,dense,0,retcode)
  else
    DAESolution(u,du,nothing,nothing,t,interp_data,prob,alg,interp,dense,0,retcode)
  end
end

function build_solution(sol::AbstractDAESolution,u_analytic,errors)
  DAESolution(sol.u,sol.du,u_analytic,errors,sol.t,sol.interp_data,
                    sol.prob,sol.alg,sol.interp,sol.dense,sol.tslocation,sol.retcode)
end
