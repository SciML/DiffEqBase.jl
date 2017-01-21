resize!(i::DEIntegrator,ii::Int) = error("This method has not been implemented for the integrator")
deleteat!(i::DEIntegrator,ii::Int) = error("This method has not been implemented for the integrator")
u_cache(i::DEIntegrator) = error("This method has not been implemented for the integrator")
du_cache(i::DEIntegrator) = error("This method has not been implemented for the integrator")
full_cache(i::DEIntegrator) = error("This method has not been implemented for the integrator")
terminate!(i::DEIntegrator) = error("This method has not been implemented for the integrator")
get_du(i::DEIntegrator) = error("This method has not been implemented for the integrator")
get_dt(i::DEIntegrator) = error("This method has not been implemented for the integrator")
get_proposed_dt(i::DEIntegrator) = error("This method has not been implemented for the integrator")
modify_proposed_dt!(i::DEIntegrator) = error("This method has not been implemented for the integrator")
u_modified!(i::DEIntegrator,bool) = error("This method has not been implemented for the integrator")
savevalues!(i::DEIntegrator) = error("This method has not been implemented for the integrator")
add_tstop!(i::DEIntegrator,t) = error("This method has not been implemented for the integrator")
add_saveat!(i::DEIntegrator,t) = error("This method has not been implemented for the integrator")
set_abstol!(i::DEIntegrator,t) = error("This method has not been implemented for the integrator")
set_reltol!(i::DEIntegrator,t) = error("This method has not been implemented for the integrator")

### Abstract Interface

immutable IntegratorTuples{I}
 integrator::I
end

start(tup::IntegratorTuples) = start(tup.integrator)

function next(tup::IntegratorTuples,state)
  state += 1
  step!(tup.integrator) # Iter updated in the step! header
  # Next is callbacks -> iterator  -> top
  (tup.integrator.t,tup.integrator.u),state
end

done(tup::IntegratorTuples,state) = done(tup.integrator,state)

tuples(integrator::DEIntegrator) = IntegratorTuples(integrator)

immutable IntegratorIntervals{I}
 integrator::I
end

start(tup::IntegratorIntervals) = start(tup.integrator)

function next(tup::IntegratorIntervals,state)
  state += 1
  step!(tup.integrator) # Iter updated in the step! header
  # Next is callbacks -> iterator  -> top
  (tup.integrator.tprev,tup.integrator.uprev,tup.integrator.t,tup.integrator.u),state
end

done(tup::IntegratorIntervals,state) = done(tup.integrator,state)

intervals(integrator::DEIntegrator) = IntegratorIntervals(integrator)

@recipe function f(integrator::DEIntegrator;
                    denseplot=integrator.opts.calck && integrator.iter>0,
                    plotdensity =10,
                    plot_analytic=false,vars=nothing)

  vars = interpret_vars(vars,integrator.sol)

  if denseplot
    # Generate the points from the plot from dense function
    plott = collect(Ranges.linspace(integrator.tprev,integrator.t,plotdensity))
    plot_timeseries = integrator(plott)
    if plot_analytic
      plot_analytic_timeseries = [integrator.sol.prob.analytic(t,integrator.sol.prob.u0) for t in plott]
    end
  end # if not denseplot, we'll just get the values right from the integrator.

  dims = length(vars[1])
  for var in vars
    @assert length(var) == dims
  end
  # Should check that all have the same dims!


  plot_vecs = []
  for i in 1:dims
    push!(plot_vecs,[])
  end

  labels = String[]# Array{String, 2}(1, length(vars)*(1+plot_analytic))
  for x in vars
    for j in 1:dims
      if denseplot
        push!(plot_vecs[j], u_n(plot_timeseries, x[j],integrator.sol,plott,plot_timeseries))
      else # just get values
        if x[j] == 0
          push!(plot_vecs[j], integrator.t)
        elseif x[j]==1 && !(typeof(integrator.u) <: AbstractArray)
          push!(plot_vecs[j], integrator.u)
        else
          push!(plot_vecs[j], integrator.u[x[j]])
        end
      end
    end
    add_labels!(labels,x,dims,integrator.sol)
  end

  if plot_analytic
    for x in vars
      for j in 1:dims
        if denseplot
          push!(plot_vecs[j], u_n(plot_timeseries, x[j],sol,plott,plot_timeseries))
        else # Just get values
          if x[j] == 0
            push!(plot_vecs[j], integrator.t)
          elseif x[j]==1 && !(typeof(integrator.u) <: AbstractArray)
            push!(plot_vecs[j], integrator.sol.prob.analytic(integrator.t,integrator.sol[1]))
          else
            push!(plot_vecs[j], integrator.sol.prob.analytic(integrator.t,integrator.sol[1])[x[j]])
          end
        end
      end
      add_labels!(labels,x,dims,integrator.sol)
    end
  end

  xflip --> integrator.tdir < 0

  if denseplot
    seriestype --> :line
  else
    seriestype --> :scatter
  end

  linewidth --> 3
  #xtickfont --> font(11)
  #ytickfont --> font(11)
  #legendfont --> font(11)
  #guidefont  --> font(11)
  label --> reshape(labels,1,length(labels))
  (plot_vecs...)
end
