@recipe function f(sol::AbstractODESolution;sensitivity=false,plot_analytic=false,denseplot=true,plotdensity=100)
  plotseries = Vector{Any}(0)
  if typeof(sol) <: AbstractSDESolution; denseplot=false; end

  if denseplot && sol.dense # Generate the points from the plot from dense function
    plott = collect(Ranges.linspace(sol.t[1],sol.t[end],plotdensity))
    plot_timeseries = sol(plott)
    if plot_analytic
      plot_analytic_timeseries = Vector{typeof(sol.u)}(length(plott))
      for i in eachindex(plott)
        tmp[i] = sol.prob.analytic(plott[i],sol.prob.u0)
      end
    end
  else # Plot for not dense output use the timeseries itself
    plot_timeseries = sol.timeseries
    if plot_analytic
      plot_analytic_timeseries = sol.timeseries_analytic
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
    if typeof(sol.u) <: AbstractArray
      for i in eachindex(sol.u)
        tmp = Vector{eltype(sol.u)}(length(plot_timeseries))
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
  lw --> 3
  #xtickfont --> font(11)
  #ytickfont --> font(11)
  #legendfont --> font(11)
  #guidefont  --> font(11)
  plott, plotseries
end
