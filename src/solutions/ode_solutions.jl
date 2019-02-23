struct ODESolution{T,N,uType,uType2,DType,tType,rateType,P,A,IType} <: AbstractODESolution{T,N}
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
  destats::DEStats
  retcode::Symbol
end
(sol::ODESolution)(t,deriv::Type=Val{0};idxs=nothing,continuity=:left) = sol.interp(t,idxs,deriv,sol.prob.p,continuity)
(sol::ODESolution)(v,t,deriv::Type=Val{0};idxs=nothing,continuity=:left) = sol.interp(v,t,idxs,deriv,sol.prob.p,continuity)

function build_solution(
        prob::Union{AbstractODEProblem,AbstractDDEProblem},
        alg,t,u;timeseries_errors=length(u)>2,
        dense=false,dense_errors=dense,
        calculate_error = true,
        k=[],
        du=[],
        interp = !isempty(du) ? HermiteInterpolation(t,u,du) : LinearInterpolation(t,u),
        retcode = :Default, destats=DEStats(), kwargs...)

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
    if !(typeof(prob.u0) <: Tuple)
      u_analytic = Vector{typeof(prob.u0)}()
      errors = Dict{Symbol,real(eltype(prob.u0))}()
    else
      u_analytic = Vector{typeof(ArrayPartition(prob.u0))}()
      errors = Dict{Symbol,real(eltype(prob.u0[1]))}()
    end

    sol = ODESolution{T,N,typeof(u),typeof(u_analytic),typeof(errors),typeof(t),typeof(k),
                      typeof(prob),typeof(alg),typeof(interp)}(u,u_analytic,
                       errors,t,k,prob,alg,interp,dense,0,destats,retcode)
    if calculate_error
      calculate_solution_errors!(sol;timeseries_errors=timeseries_errors,dense_errors=dense_errors)
    end
    return sol
  else
    return ODESolution{T,N,typeof(u),Nothing,Nothing,typeof(t),typeof(k),
                       typeof(prob),typeof(alg),typeof(interp)}(u,nothing,nothing,
                       t,k,prob,alg,interp,dense,0,destats,retcode)
  end
end

function calculate_solution_errors!(sol::AbstractODESolution;fill_uanalytic=true,timeseries_errors=true,dense_errors=true)

  f = sol.prob.f

  if fill_uanalytic
    for i in 1:size(sol.u,1)
      push!(sol.u_analytic,f.analytic(sol.prob.u0,sol.prob.p,sol.t[i]))
    end
  end

  save_everystep = length(sol.u) > 2
  if !isempty(sol.u_analytic)
    sol.errors[:final] = norm(recursive_mean(abs.(sol.u[end].-sol.u_analytic[end])))

    if save_everystep && timeseries_errors
      sol.errors[:l∞] = norm(maximum(vecvecapply((x)->abs.(x),sol.u-sol.u_analytic)))
      sol.errors[:l2] = norm(sqrt(recursive_mean(vecvecapply((x)->float.(x).^2,sol.u-sol.u_analytic))))
      if sol.dense && dense_errors
        densetimes = collect(range(sol.t[1], stop=sol.t[end], length=100))
        interp_u = sol(densetimes)
        interp_analytic = VectorOfArray([f.analytic(sol.prob.u0,sol.prob.p,t) for t in densetimes])
        sol.errors[:L∞] = norm(maximum(abs.(interp_u.-interp_analytic)))
        sol.errors[:L2] = norm(sqrt(mean((interp_u.-interp_analytic).^2)))
      end
    end
  end
end

function build_solution(sol::AbstractODESolution{T,N},u_analytic,errors) where {T,N}
  ODESolution{T,N,typeof(sol.u),typeof(u_analytic),typeof(errors),typeof(sol.t),typeof(sol.k),
                     typeof(sol.prob),typeof(sol.alg),typeof(sol.interp)}(
                     sol.u,u_analytic,errors,sol.t,sol.k,sol.prob,
                     sol.alg,sol.interp,sol.dense,sol.tslocation,sol.destats,sol.retcode)
end

function solution_new_retcode(sol::AbstractODESolution{T,N},retcode) where {T,N}
  ODESolution{T,N,typeof(sol.u),typeof(sol.u_analytic),typeof(sol.errors),
                     typeof(sol.t),typeof(sol.k),
                     typeof(sol.prob),typeof(sol.alg),typeof(sol.interp)}(
                     sol.u,sol.u_analytic,sol.errors,sol.t,sol.k,sol.prob,
                     sol.alg,sol.interp,sol.dense,sol.tslocation,sol.destats,retcode)
 end

 function solution_new_tslocation(sol::AbstractODESolution{T,N},tslocation) where {T,N}
   ODESolution{T,N,typeof(sol.u),typeof(sol.u_analytic),typeof(sol.errors),
                      typeof(sol.t),typeof(sol.k),
                      typeof(sol.prob),typeof(sol.alg),typeof(sol.interp)}(
                      sol.u,sol.u_analytic,sol.errors,sol.t,sol.k,sol.prob,
                      sol.alg,sol.interp,sol.dense,tslocation,sol.destats,sol.retcode)
  end

  function solution_slice(sol::AbstractODESolution{T,N},I) where {T,N}
    ODESolution{T,N,typeof(sol.u),typeof(sol.u_analytic),typeof(sol.errors),
                       typeof(sol.t),typeof(sol.k),
                       typeof(sol.prob),typeof(sol.alg),typeof(sol.interp)}(
                       sol.u[I],
                       sol.u_analytic === nothing ? nothing : sol.u_analytic[I],
                       sol.errors,sol.t[I],sol.k[I],sol.prob,
                       sol.alg,sol.interp,false,sol.tslocation,sol.destats,sol.retcode)
   end
