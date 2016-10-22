Base.length(sol::DESolution) = length(sol.t)
Base.size(sol::DESolution) = (length(sol.t),size(sol.u))
Base.endof(sol::DESolution) = length(sol)
Base.getindex(sol::DESolution,i::Int) = sol.timeseries[i]
Base.getindex(sol::DESolution,i::Int,I::Int...) = sol.timeseries[i][I...]
Base.getindex(sol::DESolution,::Colon) = sol.timeseries

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
  sol.timeseries!=[] && println(io,"timeseries: $(sol.timeseries)")
  sol.trueknown && sol.timeseries_analytic!=[] && println(io,"timeseries_analytic: $(sol.timeseries_analytic)")
  nothing
end

function show(io::IO,sol::DESolution)
  print(io,"$(typeof(sol)), $(length(sol)) timesteps, final value $(sol.u)")
end
