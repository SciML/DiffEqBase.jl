### Concrete Types

type RODESolution{T,N,uType,uType2,DType,tType,randType,P,A,IType} <: AbstractRODESolution{T,N}
  u::uType
  u_analytic::uType2
  errors::DType
  t::tType
  W::randType
  prob::P
  alg::A
  interp::IType
  dense::Bool
  tslocation::Int
  retcode::Symbol
  seed::UInt64
end
(sol::RODESolution)(t,deriv::Type=Val{0};idxs=nothing) = sol.interp(t,idxs,deriv)
(sol::RODESolution)(v,t,deriv::Type=Val{0};idxs=nothing) = sol.interp(v,t,idxs,deriv)

function build_solution(
        prob::AbstractRODEProblem,
        alg,t,u;W=[],timeseries_errors = length(u) > 2,
        dense = false,dense_errors=dense,calculate_error=true,
        interp = LinearInterpolation(t,u),
        retcode = :Default,
        seed = UInt64(0), kwargs...)

  T = eltype(eltype(u))
  if typeof(prob.u0) <: Tuple
    N = length((size(ArrayPartition(prob.u0))..., length(u)))
  else
    N = length((size(prob.u0)..., length(u)))
  end

  if typeof(prob.f) <: Tuple
    f = prob.f[1]
  else
    f = prob.f
  end

  if has_analytic(f)
    u_analytic = Vector{typeof(prob.u0)}(0)
    errors = Dict{Symbol,eltype(prob.u0)}()
    sol = RODESolution{T,N,typeof(u),typeof(u_analytic),typeof(errors),typeof(t),typeof(W),
                       typeof(prob),typeof(alg),typeof(interp)}(
                       u,u_analytic,errors,t,W,prob,alg,interp,dense,0,retcode,seed)

    if calculate_error
      calculate_solution_errors!(sol;timeseries_errors=timeseries_errors,dense_errors=dense_errors)
    end

    return sol
  else
    return RODESolution{T,N,typeof(u),Void,Void,typeof(t),
                        typeof(W),typeof(prob),typeof(alg),typeof(interp)}(
                        u,nothing,nothing,t,W,prob,alg,interp,dense,0,retcode,seed)
  end
end

function calculate_solution_errors!(sol::AbstractRODESolution;fill_uanalytic=true,
                                    timeseries_errors=true,dense_errors=true)

  if typeof(sol.prob.f) <: Tuple
    f = sol.prob.f[1]
  else
    f = sol.prob.f
  end

  if fill_uanalytic
    for i in 1:length(sol)
      push!(sol.u_analytic,f(Val{:analytic},sol.t[i],sol.prob.u0,sol.W[i]))
    end
  end

  if !isempty(sol.u_analytic)
    sol.errors[:final] = recursive_mean(abs.(sol.u[end]-sol.u_analytic[end]))
    if timeseries_errors
      sol.errors[:l∞] = maximum(vecvecapply((x)->abs.(x),sol.u-sol.u_analytic))
      sol.errors[:l2] = sqrt(recursive_mean(vecvecapply((x)->float.(x).^2,sol.u-sol.u_analytic)))
    end
    if dense_errors
      densetimes = collect(linspace(sol.t[1],sol.t[end],100))
      interp_u = sol(densetimes)
      interp_analytic = [f(Val{:analytic},t,sol.u[1],sol.W(t)[1]) for t in densetimes]
      sol.errors[:L∞] = maximum(vecvecapply((x)->abs.(x),interp_u-interp_analytic))
      sol.errors[:L2] = sqrt(recursive_mean(vecvecapply((x)->float.(x).^2,interp_u-interp_analytic)))
    end
  end
end

function build_solution(sol::AbstractRODESolution,u_analytic,errors)
  T = eltype(eltype(sol.u))

  if typeof(sol.u) <: Tuple
    N = length((size(ArrayPartition(sol.u))..., length(sol.u)))
  else
    N = length((size(sol.u[1])..., length(sol.u)))
  end

  RODESolution{T,N,typeof(sol.u),typeof(u_analytic),typeof(errors),typeof(sol.t),
               typeof(sol.W),typeof(sol.prob),typeof(sol.alg),typeof(sol.interp)}(
               sol.u,u_analytic,errors,sol.t,sol.W,sol.prob,sol.alg,sol.interp,
               sol.dense,sol.tslocation,sol.retcode,sol.seed)
end
