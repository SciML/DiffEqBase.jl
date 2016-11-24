### Abstract Interface

Base.length(sol::DESolution) = length(sol.u) # Must be on u for the test solutions!
Base.endof(sol::DESolution) = length(sol)
Base.getindex(sol::DESolution,i::Int) = sol.u[i]
Base.getindex(sol::DESolution,i::Int,I::Int...) = sol.u[i][I...]
Base.getindex(sol::DESolution,::Colon) = sol.u
Base.getindex(sol::DESolution,::Colon,i::Int...) = [sol.u[j][i...] for j in eachindex(sol)]
eachindex(sol::DESolution) = eachindex(sol.t)
tuples(sol::DESolution) = tuple.(sol.t,sol.u)

function start(sol::DESolution)
  sol.tslocation = 1
  1
end

function next(sol::DESolution,state)
  state += 1
  sol.tslocation = state
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

#=
function print(io::IO, sol::DESolution)
  println(io,"$(typeof(sol))")
  println(io,"u: $(sol.u)")
  println(io,"t: $(sol.t)")
  nothing
end
=#

#=
function show(io::IO,sol::DESolution)
  print(io,"$(typeof(sol))")
end
=#

@recipe function f(sol::AbstractODESolution;plot_analytic=false,denseplot=true,plotdensity=100)
  plotseries = Vector{Any}(0)
  if typeof(sol) <: AbstractSDESolution; denseplot=false; end

  if denseplot && sol.dense # Generate the points from the plot from dense function
    plott = collect(Ranges.linspace(sol.t[1],sol.t[end],plotdensity))
    plot_timeseries = sol(plott)
    if plot_analytic
      plot_analytic_timeseries = [sol.prob.analytic(t,sol.prob.u0) for t in plott]
    end
  else # Plot for not dense output use the timeseries itself
    plot_timeseries = sol.u
    if plot_analytic
      plot_analytic_timeseries = sol.u_analytic
    end
    plott = sol.t
  end

  # Make component-wise plots
  if typeof(sol[1]) <:AbstractArray
    for i in eachindex(sol[1])
      tmp = Vector{eltype(sol[1])}(length(plot_timeseries))
      for j in 1:length(plot_timeseries)
        tmp[j] = plot_timeseries[j][i]
      end
      push!(plotseries,tmp)
    end
  else
    push!(plotseries,plot_timeseries)
  end
  if plot_analytic
    if typeof(sol[1]) <: AbstractArray
      for i in eachindex(sol[1])
        tmp = Vector{eltype(sol[1])}(length(plot_timeseries))
        for j in 1:length(plot_timeseries)
          tmp[j] = plot_analytic_timeseries[j][i]
        end
        push!(plotseries,tmp)
      end
    else
      push!(plotseries,plot_analytic_timeseries)
    end
  end

  seriestype --> :path
  linewidth --> 3
  #xtickfont --> font(11)
  #ytickfont --> font(11)
  #legendfont --> font(11)
  #guidefont  --> font(11)
  plott, plotseries
end
