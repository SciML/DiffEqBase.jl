resize!(i::DEIntegrator,ii::Int) = error("resize!: method has not been implemented for the integrator")
deleteat!(i::DEIntegrator,ii) = error("deleteat!: method has not been implemented for the integrator")
addat!(i::DEIntegrator,ii,val=zeros(length(idxs))) = error("addat!: method has not been implemented for the integrator")
user_cache(i::DEIntegrator) = error("user_cache: method has not been implemented for the integrator")
u_cache(i::DEIntegrator) = error("u_cache: method has not been implemented for the integrator")
du_cache(i::DEIntegrator) = error("du_cache: method has not been implemented for the integrator")
full_cache(i::DEIntegrator) = error("full_cache: method has not been implemented for the integrator")
resize_non_user_cache!(i::DEIntegrator,ii::Int) = error("resize_non_user_cache!: method has not been implemented for the integrator")
deleteat_non_user_cache!(i::DEIntegrator,idxs) = error("deleteat_non_user_cache!: method has not been implemented for the integrator")
addat_non_user_cache!(i::DEIntegrator,idxs) = error("addat_non_user_cache!: method has not been implemented for the integrator")
terminate!(i::DEIntegrator) = error("terminate!: method has not been implemented for the integrator")
get_du(i::DEIntegrator) = error("get_du: method has not been implemented for the integrator")
get_dt(i::DEIntegrator) = error("get_dt: method has not been implemented for the integrator")
get_proposed_dt(i::DEIntegrator) = error("get_proposed_dt: method has not been implemented for the integrator")
modify_proposed_dt!(i::DEIntegrator) = error("modify_proposed_dt!: method has not been implemented for the integrator")
u_modified!(i::DEIntegrator,bool) = error("u_modified!: method has not been implemented for the integrator")
savevalues!(i::DEIntegrator) = error("savevalues!: method has not been implemented for the integrator")
add_tstop!(i::DEIntegrator,t) = error("add_tstop!: method has not been implemented for the integrator")
add_saveat!(i::DEIntegrator,t) = error("add_saveat!: method has not been implemented for the integrator")
set_abstol!(i::DEIntegrator,t) = error("set_abstol!: method has not been implemented for the integrator")
set_reltol!(i::DEIntegrator,t) = error("set_reltol!: method has not been implemented for the integrator")

### Addat isn't a real thing. Let's make it a real thing Gretchen

function addat!(a::AbstractArray,idxs,val=zeros(length(idxs)))
  flip_range = UnitRange(idxs.start,idxs.start-length(idxs))
  splice!(a,flip_range,val)
end

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
                    denseplot=(integrator.opts.calck || typeof(integrator) <: AbstractSDEIntegrator)  && integrator.iter>0,
                    plotdensity =10,
                    plot_analytic=false,vars=nothing)

  int_vars = interpret_vars(vars,integrator.sol)

  if denseplot
    # Generate the points from the plot from dense function
    plott = collect(Ranges.linspace(integrator.tprev,integrator.t,plotdensity))
    plot_timeseries = integrator(plott)
    if plot_analytic
      plot_analytic_timeseries = [integrator.sol.prob.analytic(t,integrator.sol.prob.u0) for t in plott]
    end
  end # if not denseplot, we'll just get the values right from the integrator.

  dims = length(int_vars[1])
  for var in int_vars
    @assert length(var) == dims
  end
  # Should check that all have the same dims!


  plot_vecs = []
  for i in 1:dims
    push!(plot_vecs,[])
  end

  labels = String[]# Array{String, 2}(1, length(int_vars)*(1+plot_analytic))
  for x in int_vars
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
    for x in int_vars
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
    seriestype --> :path
  else
    seriestype --> :scatter
  end

  if typeof(vars) <: Tuple && eltype(vars) == Symbol
    xlabel --> vars[1]
    ylabel --> vars[2]
    if length(vars) > 2
      zlabel --> vars[3]
    end
  end
  if first.(int_vars) == zeros(length(int_vars))
    xlabel --> "t"
  end

  linewidth --> 3
  #xtickfont --> font(11)
  #ytickfont --> font(11)
  #legendfont --> font(11)
  #guidefont  --> font(11)
  label --> reshape(labels,1,length(labels))
  (plot_vecs...)
end
