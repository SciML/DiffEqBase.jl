struct DAESolution{T,N,uType,duType,uType2,DType,tType,P,A,ID} <: AbstractDAESolution{T,N}
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
(sol::DAESolution)(t,deriv::Type=Val{0};idxs=nothing) = sol.interp(t,idxs,deriv,sol.prob.p)
(sol::DAESolution)(v,t,deriv::Type=Val{0};idxs=nothing) = sol.interp(v,t,idxs,deriv,sol.prob.p)

function build_solution(
        prob::AbstractDAEProblem{uType,duType,tType,isinplace},
alg,t,u;dense=false,du=[],
interp = !isempty(du) ? HermiteInterpolation(t,u,du) : LinearInterpolation(t,u),
timeseries_errors=true,dense_errors=true, retcode = :Default, kwargs...) where {uType,duType,tType,isinplace}

  T = eltype(eltype(u))
  N = length((size(u[1])..., length(u)))

  if has_analytic(prob.f)
    u_analytic = Vector{typeof(prob.u0)}(0)
    for i in 1:size(u,1)
      push!(u_analytic,prob.analytic(prob.u0,t[i]))
    end

    save_everystep = length(u) > 2

    errors = Dict{Symbol,real(eltype(prob.u0))}()
    if !isempty(u_analytic)
      errors[:final] = norm(recursive_mean(abs.(u[end]-u_analytic[end])))

      if save_everystep && timeseries_errors
        errors[:l∞] = norm(maximum(vecvecapply((x)->abs.(x),u-u_analytic)))
        errors[:l2] = norm(sqrt(recursive_mean(vecvecapply((x)->float.(x).^2,u-u_analytic))))
        if dense && dense_errors
          densetimes = collect(range(t[1], stop=t[end], length=100))
          interp_u = interp(densetimes)
          interp_analytic = [prob.analytic(t,u[1]) for t in densetimes]
          errors[:L∞] = norm(maximum(vecvecapply((x)->abs.(x),interp_u-interp_analytic)))
          errors[:L2] = norm(sqrt(recursive_mean(vecvecapply((x)->float.(x).^2,interp_u-interp_analytic))))
        end
      end
    end
    DAESolution{T,N,typeof(u),typeof(du),typeof(u_analytic),typeof(errors),typeof(t),
                       typeof(prob),typeof(alg),typeof(interp)}(u,du,u_analytic,errors,t,prob,alg,interp,dense,0,retcode)
  else
    DAESolution{T,N,typeof(u),typeof(du),Nothing,Nothing,typeof(t),
                       typeof(prob),typeof(alg),typeof(interp)}(u,du,nothing,nothing,t,prob,alg,interp,dense,0,retcode)
  end
end

function build_solution(sol::AbstractDAESolution{T,N},u_analytic,errors) where {T,N}
  DAESolution{T,N,typeof(sol.u),typeof(sol.du),typeof(u_analytic),typeof(errors),typeof(sol.t),
                     typeof(sol.prob),typeof(sol.alg),typeof(sol.interp)}(
                     sol.u,sol.du,u_analytic,errors,sol.t,
              sol.prob,sol.alg,sol.interp,sol.dense,sol.tslocation,sol.retcode)
end

function solution_new_retcode(sol::AbstractDAESolution{T,N},retcode) where {T,N}
  DAESolution{T,N,typeof(sol.u),typeof(sol.du),typeof(sol.u_analytic),
              typeof(sol.errors),typeof(sol.t),
              typeof(sol.prob),typeof(sol.alg),typeof(sol.interp)}(
              sol.u,sol.du,sol.u_analytic,sol.errors,sol.t,
              sol.prob,sol.alg,sol.interp,sol.dense,sol.tslocation,retcode)
end

function solution_new_tslocation(sol::AbstractDAESolution{T,N},tslocation) where {T,N}
  DAESolution{T,N,typeof(sol.u),typeof(sol.du),typeof(sol.u_analytic),
              typeof(sol.errors),typeof(sol.t),
              typeof(sol.prob),typeof(sol.alg),typeof(sol.interp)}(
              sol.u,sol.du,sol.u_analytic,sol.errors,sol.t,
              sol.prob,sol.alg,sol.interp,sol.dense,tslocation,sol.retcode)
end
