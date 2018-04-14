"""
    step!(integ::DEIntegrator [, dt [, stop_at_tdt]])
Perform one (successful) step on the integrator.

Alternative, if a `dt` is given, then `step!` the integrator until
there is a temporal difference `≥ dt` in `integ.t`.  It returns the
actual temporal difference advanced by the integrator.  When `true` is
passed to the optional third argument, the integrator advances exactly
`dt`.
"""
function step!(d::DEIntegrator) error("Integrator stepping is not implemented") end
resize!(i::DEIntegrator,ii::Int) = error("resize!: method has not been implemented for the integrator")
deleteat!(i::DEIntegrator,ii) = error("deleteat!: method has not been implemented for the integrator")
addat!(i::DEIntegrator,ii,val=zeros(length(idxs))) = error("addat!: method has not been implemented for the integrator")
get_tmp_cache(i::DEIntegrator) = error("get_tmp_cache!: method has not been implemented for the integrator")
user_cache(i::DEIntegrator) = error("user_cache: method has not been implemented for the integrator")
u_cache(i::DEIntegrator) = error("u_cache: method has not been implemented for the integrator")
du_cache(i::DEIntegrator) = error("du_cache: method has not been implemented for the integrator")
full_cache(i::DEIntegrator) = error("full_cache: method has not been implemented for the integrator")
resize_non_user_cache!(i::DEIntegrator,ii::Int) = error("resize_non_user_cache!: method has not been implemented for the integrator")
deleteat_non_user_cache!(i::DEIntegrator,idxs) = error("deleteat_non_user_cache!: method has not been implemented for the integrator")
addat_non_user_cache!(i::DEIntegrator,idxs) = error("addat_non_user_cache!: method has not been implemented for the integrator")
terminate!(i::DEIntegrator) = error("terminate!: method has not been implemented for the integrator")
get_du(i::DEIntegrator) = error("get_du: method has not been implemented for the integrator")
get_du!(out,i::DEIntegrator) = error("get_du: method has not been implemented for the integrator")
get_dt(i::DEIntegrator) = error("get_dt: method has not been implemented for the integrator")
get_proposed_dt(i::DEIntegrator) = error("get_proposed_dt: method has not been implemented for the integrator")
set_proposed_dt!(i::DEIntegrator) = error("modify_proposed_dt!: method has not been implemented for the integrator")
u_modified!(i::DEIntegrator,bool) = error("u_modified!: method has not been implemented for the integrator")
savevalues!(i::DEIntegrator) = error("savevalues!: method has not been implemented for the integrator")
add_tstop!(i::DEIntegrator,t) = error("add_tstop!: method has not been implemented for the integrator")
add_saveat!(i::DEIntegrator,t) = error("add_saveat!: method has not been implemented for the integrator")
set_abstol!(i::DEIntegrator,t) = error("set_abstol!: method has not been implemented for the integrator")
set_reltol!(i::DEIntegrator,t) = error("set_reltol!: method has not been implemented for the integrator")
reinit!(integrator::DEIntegrator,args...; kwargs...) =
       error("reinit!: method has not been implemented for the integrator")
auto_dt_reset!(integrator::DEIntegrator) = error("auto_dt_reset!: method has not been implemented for the integrator")

"""
    set_t!(integrator::DEIntegrator, t::Real)

Set current time point of the `integrator` to `t`.
"""
set_t!(integrator::DEIntegrator, t::Real) =
    error("set_t!: method has not been implemented for the integrator")

"""
    set_u!(integrator::DEIntegrator, u)

Set current state of the `integrator` to `u`.
"""
set_u!(integrator::DEIntegrator, u) =
    error("set_u!: method has not been implemented for the integrator")


"""
    set_ut!(integrator::DEIntegrator, u, t)

Set current state of the `integrator` to `u` and `t`
"""
function set_ut!(integrator::DEIntegrator, u, t)
  DiffEqBase.set_u!(integrator, u)
  DiffEqBase.set_t!(integrator, t)
end

### Addat isn't a real thing. Let's make it a real thing Gretchen

function addat!(a::AbstractArray,idxs,val=zeros(length(idxs)))
  flip_range = UnitRange(idxs.start,idxs.start-length(idxs))
  splice!(a,flip_range,val)
end

### Integrator traits

has_reinit(i::DEIntegrator) = false

### Display

Base.summary(I::DEIntegrator) = string("Integrator with uType ",typeof(I.u)," and tType ",typeof(I.t))
function Base.show(io::IO, A::DEIntegrator)
  println(io,string("t: ",A.t))
  print(io,"u: ")
  show(io, A.u)
end
function Base.show(io::IO, m::MIME"text/plain", A::DEIntegrator)
  println(io,string("t: ",A.t))
  print(io,"u: ")
  show(io,m,A.u)
end

### Error check (retcode)

last_step_failed(integrator::DEIntegrator) = false

"""
    check_error(integrator)

Check state of `integrator` and return one of the
[Return Codes](http://docs.juliadiffeq.org/latest/basics/solution.html#Return-Codes-(RetCodes)-1)
"""
function check_error(integrator::DEIntegrator)
  # This implementation is intended to be used for ODEIntegrator and
  # SDEIntegrator.
  if integrator.iter > integrator.opts.maxiters
    if integrator.opts.verbose
      warn("Interrupted. Larger maxiters is needed.")
    end
    return :MaxIters
  end
  if !integrator.opts.force_dtmin && integrator.opts.adaptive && abs(integrator.dt) <= abs(integrator.opts.dtmin)
    if integrator.opts.verbose
      warn("dt <= dtmin. Aborting. If you would like to force continuation with dt=dtmin, set force_dtmin=true")
    end
    return :DtLessThanMin
  end
  if integrator.opts.unstable_check(integrator.dt,integrator.u,integrator.p,integrator.t)
    if integrator.opts.verbose
      warn("Instability detected. Aborting")
    end
    return :Unstable
  end
  if last_step_failed(integrator)
    if integrator.opts.verbose
      warn("Newton steps could not converge and algorithm is not adaptive. Use a lower dt.")
    end
    return :ConvergenceFailure
  end
  return :Success
end

function postamble! end

"""
    check_error!(integrator)

Same as `check_error` but also set solution's return code
(`integrator.sol.retcode`) and run `postamble!`.
"""
function check_error!(integrator::DEIntegrator)
  code = check_error(integrator)
  if code != :Success
    integrator.sol = solution_new_retcode(integrator.sol, code)
    postamble!(integrator)
  end
  return code
end

### Default Iterator Interface

function start(integrator::DEIntegrator)
  0
end

@inline function next(integrator::DEIntegrator,state)
  state += 1
  step!(integrator) # Iter updated in the step! header
  # Next is callbacks -> iterator  -> top
  integrator,state
end

@inline function done(integrator::DEIntegrator, _)
  if check_error!(integrator) != :Success
    return true
  elseif isempty(integrator.opts.tstops)
    postamble!(integrator)
    return true
  elseif integrator.just_hit_tstop
    integrator.just_hit_tstop = false
    if integrator.opts.stop_at_next_tstop
      postamble!(integrator)
      return true
    end
  end
  false
end

done(integrator::DEIntegrator) = done(integrator,integrator.iter)

eltype(integrator::DEIntegrator) = typeof(integrator)

### Abstract Interface

struct IntegratorTuples{I}
 integrator::I
end

start(tup::IntegratorTuples) = start(tup.integrator)

function next(tup::IntegratorTuples,state)
  state += 1
  step!(tup.integrator) # Iter updated in the step! header
  # Next is callbacks -> iterator  -> top
  (tup.integrator.u,tup.integrator.t),state
end

done(tup::IntegratorTuples,state) = done(tup.integrator,state)

tuples(integrator::DEIntegrator) = IntegratorTuples(integrator)

struct IntegratorIntervals{I}
 integrator::I
end

start(tup::IntegratorIntervals) = start(tup.integrator)

function next(tup::IntegratorIntervals,state)
  state += 1
  step!(tup.integrator) # Iter updated in the step! header
  # Next is callbacks -> iterator  -> top
  (tup.integrator.uprev,tup.integrator.tprev,tup.integrator.u,tup.integrator.t),state
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
    plott = collect(linspace(integrator.tprev,integrator.t,plotdensity))
    plot_timeseries = integrator(plott)
    if plot_analytic
      plot_analytic_timeseries = [integrator.sol.prob.f(Val{:analytic},t,integrator.sol.prob.u0) for t in plott]
    end
  end # if not denseplot, we'll just get the values right from the integrator.

  dims = length(int_vars[1])
  for var in int_vars
    @assert length(var) == dims
  end
  # Should check that all have the same dims!


  plot_vecs = []
  for i in 2:dims
    push!(plot_vecs,[])
  end

  labels = String[]# Array{String, 2}(1, length(int_vars)*(1+plot_analytic))
  for x in int_vars
    for j in 2:dims
      if denseplot
        push!(plot_vecs[j-1], u_n(plot_timeseries, x[j],integrator.sol,plott,plot_timeseries))
      else # just get values
        if x[j] == 0
          push!(plot_vecs[j-1], integrator.t)
        elseif x[j]==1 && !(typeof(integrator.u) <: AbstractArray)
          push!(plot_vecs[j-1], integrator.u)
        else
          push!(plot_vecs[j-1], integrator.u[x[j]])
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
            push!(plot_vecs[j], integrator.sol.prob.f(Val{:analytic},integrator.t,integrator.sol[1]))
          else
            push!(plot_vecs[j], integrator.sol.prob.f(Val{:analytic},integrator.t,integrator.sol[1])[x[j]])
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

  # Special case labels when vars = (:x,:y,:z) or (:x) or [:x,:y] ...
  if typeof(vars) <: Tuple && (typeof(vars[1]) == Symbol && typeof(vars[2]) == Symbol)
    xlabel --> vars[1]
    ylabel --> vars[2]
    if length(vars) > 2
      zlabel --> vars[3]
    end
  end
  if getindex.(int_vars,1) == zeros(length(int_vars)) || getindex.(int_vars,2) == zeros(length(int_vars))
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

function step!(integ::DEIntegrator, dt::Real, stop_at_tdt = false)
    (dt * integ.tdir) < 0 && error("Cannot step backward.")
    t = integ.t
    next_t = t+dt
    stop_at_tdt && add_tstop!(integ,next_t)
    while integ.t < next_t
        step!(integ)
        integ.sol.retcode in (:Default, :Success) || break
    end
    return integ.t - t
end
