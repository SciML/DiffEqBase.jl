Base.length(sol::DESolution) = length(sol.u)
Base.endof(sol::DESolution) = length(sol)
Base.getindex(sol::DESolution,i::Int) = sol.u[i]
Base.getindex(sol::DESolution,i::Int,I::Int...) = sol.u[i][I...]
Base.getindex(sol::DESolution,::Colon) = sol.u

function start(sol::DESolution)
  #sol.tslocation = state
  1
end

function next(sol::DESolution,state)
  state += 1
  #sol.tslocation = state
  (sol,state)
end

function done(sol::DESolution,state)
  state >= length(sol)
end

function eltype(sol::DESolution)
  if typeof(sol[1]) <: AbstractArray
    return typeof(sol[1][1])
  else
    return typeof(sol[1])
  end
end

function print(io::IO, sol::DESolution)
  println(io,"$(typeof(sol)) with $(length(sol)) timesteps.")
  println(io,"u: $(sol.u)")
  println(io,"t: $(sol.t)")
  nothing
end

function show(io::IO,sol::DESolution)
  print(io,"$(typeof(sol)), $(length(sol)) timesteps, final value $(sol[end])")
end


"""
`ODESolution`

Holds the data for the solution to an ODE problem.

### Fields

* `u::Array{Float64}`: The solution
* `u_analytic::AbstractArrayOrVoid`: The true solution at the final timepoint.
* `errors`: A dictionary of the error calculations.
* `t::AbstractArrayOrVoid`: All the t's in the solution. Only saved if `save_timeseries=true`
  is specified in the solver.
* `prob::DEProblem`: Holds the problem object used to define the problem.

"""
type ODESolution{uType,uEltype,tType,rateType,P,A} <: AbstractODESolution
  u::uType
  u_analytic
  errors::Dict{Symbol,uEltype}
  t::tType
  k::rateType
  prob::P
  alg::A
  interp::Function
  dense::Bool
end

function ODESolution{uType,tType,isinplace}(t,u,
        prob::AbstractODEProblem{uType,tType,Val{isinplace}},
        alg;u_analytic=[],k=[],saveat=[],
        interp = (tvals) -> nothing,
        timeseries_errors=true,dense_errors=true)

  save_timeseries = length(u) > 2

  dense = length(k)>1
  errors = Dict{Symbol,eltype(u[1])}()
  if !isempty(u_analytic)
    errors[:final] = mean(abs.(u[end]-u_analytic[end]))

    if save_timeseries && timeseries_errors
      errors[:l∞] = maximum(vecvecapply((x)->abs.(x),u-u_analytic))
      errors[:l2] = sqrt(mean(vecvecapply((x)->float(x).^2,u-u_analytic)))
      if dense && dense_errors
        densetimes = collect(linspace(t[1],t[end],100))
        interp_u = interp(densetimes)
        interp_analytic = [prob.analytic(t,u[1]) for t in densetimes]
        errors[:L∞] = maximum(vecvecapply((x)->abs.(x),interp_u-interp_analytic))
        errors[:L2] = sqrt(mean(vecvecapply((x)->float(x).^2,interp_u-interp_analytic)))
      end
    end
  end
  return(ODESolution(u,u_analytic,errors,t,k,prob,alg,interp,dense))
end

(sol::ODESolution)(t) = sol.interp(t)
