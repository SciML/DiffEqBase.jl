type DAESolution{T,N,uType,duType,uType2,DType,tType,P,A,ID} <: AbstractODESolution{T,N}
  u::uType
  du::duType
  u_analytic::uType2
  errors::DType
  t::tType
  prob::P
  alg::A
  interp::ID
  dense::Bool
  tslocation::Int
  retcode::Symbol
end
(sol::DAESolution)(t,deriv::Type=Val{0};idxs=nothing) = sol.interp(t,idxs,deriv)
(sol::DAESolution)(v,t,deriv::Type=Val{0};idxs=nothing) = sol.interp(v,t,idxs,deriv)

function build_solution{uType,duType,tType,isinplace}(
        prob::AbstractDAEProblem{uType,duType,tType,isinplace},
        alg,t,u;dense=false,du=[],
        interp = (t,idxs,deriv) -> error("Post-solution interpolation is not compatible with this algorithm. For other choices, see the solver compatibility chart: http://docs.juliadiffeq.org/latest/basics/compatibility_chart.html"),
        timeseries_errors=true,dense_errors=true, retcode = :Default, kwargs...)

  T = eltype(eltype(u))
  N = length((size(u[1])..., length(u)))

  if has_analytic(prob.f)
    u_analytic = Vector{typeof(prob.u0)}(0)
    for i in 1:size(u,1)
      push!(u_analytic,prob.analytic(t[i],prob.u0))
    end

    save_everystep = length(u) > 2

    errors = Dict{Symbol,eltype(u[1])}()
    if !isempty(u_analytic)
      errors[:final] = recursive_mean(abs.(u[end]-u_analytic[end]))

      if save_everystep && timeseries_errors
        errors[:l∞] = maximum(vecvecapply((x)->abs.(x),u-u_analytic))
        errors[:l2] = sqrt(recursive_mean(vecvecapply((x)->float.(x).^2,u-u_analytic)))
        if dense && dense_errors
          densetimes = collect(linspace(t[1],t[end],100))
          interp_u = interp(densetimes)
          interp_analytic = [prob.analytic(t,u[1]) for t in densetimes]
          errors[:L∞] = maximum(vecvecapply((x)->abs.(x),interp_u-interp_analytic))
          errors[:L2] = sqrt(recursive_mean(vecvecapply((x)->float.(x).^2,interp_u-interp_analytic)))
        end
      end
    end
    DAESolution{T,N,typeof(u),typeof(du),typeof(u_analytic),typeof(errors),typeof(t),
                       typeof(prob),typeof(alg),typeof(interp)}(u,du,u_analytic,errors,t,prob,alg,interp,dense,0,retcode)
  else
    DAESolution{T,N,typeof(u),typeof(du),Void,Void,typeof(t),
                       typeof(prob),typeof(alg),typeof(interp)}(u,du,nothing,nothing,t,prob,alg,interp,dense,0,retcode)
  end
end

function build_solution(sol::AbstractDAESolution,u_analytic,errors)
  DAESolution{T,N,typeof(sol.u),typeof(sol.du),typeof(u_analytic),typeof(errors),typeof(sol.t),
                     typeof(sol.prob),typeof(sol.alg),typeof(sol.interp)}(
                     sol.u,sol.du,u_analytic,errors,sol.t,sol.interp_data,
              sol.prob,sol.alg,sol.interp,sol.dense,sol.tslocation,sol.retcode)
end
