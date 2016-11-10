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
  if sol.trueknown
    str="Analytical solution is known."
  else
    str="No analytical solution is known."
  end
  println(io,"$(typeof(sol)) with $(length(sol)) timesteps. $str")
  println(io,"u: $(sol.u)")
  sol.trueknown && println(io,"errors: $(sol.errors)")
  sol.t!=[] && println(io,"t: $(sol.t)")
  sol.trueknown && sol.timeseries_analytic!=[] && println(io,"timeseries_analytic: $(sol.timeseries_analytic)")
  nothing
end

function show(io::IO,sol::DESolution)
  print(io,"$(typeof(sol)), $(length(sol)) timesteps, final value $(sol[end])")
end
