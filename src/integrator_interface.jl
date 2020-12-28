"""
    step!(integ::DEIntegrator [, dt [, stop_at_tdt]])
Perform one (successful) step on the integrator.

Alternative, if a `dt` is given, then `step!` the integrator until
there is a temporal difference `≥ dt` in `integ.t`.  When `true` is
passed to the optional third argument, the integrator advances exactly
`dt`.
"""
function step!(d::DEIntegrator) error("Integrator stepping is not implemented") end

"""
    resize(integrator::DEIntegrator,k::Int)

Resizes the DE to a size `k`. This chops off the end of the array, or adds blank values at the end, depending on whether
`k > length(integrator.u)`.
"""
Base.resize!(i::DEIntegrator,ii::Int) = error("resize!: method has not been implemented for the integrator")

"""
    deleteat!(integrator::DEIntegrator,idxs)

Shrinks the ODE by deleting the `idxs` components.
"""
Base.deleteat!(i::DEIntegrator,ii) = error("deleteat!: method has not been implemented for the integrator")

"""
    addat!(integrator::DEIntegrator,idxs,val)

Grows the ODE by adding the `idxs` components. Must be contiguous indices.
"""
addat!(i::DEIntegrator,idxs,val=zeros(length(idxs))) = error("addat!: method has not been implemented for the integrator")

"""
    get_tmp_cache(i::DEIntegrator)

Returns a tuple of internal cache vectors which are safe to use as temporary arrays. This should be used
for integrator interface and callbacks which need arrays to write into in order to be non-allocating.
The length of the tuple is dependent on the method.
"""
get_tmp_cache(i::DEIntegrator) = error("get_tmp_cache!: method has not been implemented for the integrator")
user_cache(i::DEIntegrator) = error("user_cache: method has not been implemented for the integrator")
u_cache(i::DEIntegrator) = error("u_cache: method has not been implemented for the integrator")
du_cache(i::DEIntegrator) = error("du_cache: method has not been implemented for the integrator")
ratenoise_cache(i::DEIntegrator) = ()
rand_cache(i::DEIntegrator) = ()

"""
    full_cache(i::DEIntegrator)

Returns an iterator over the cache arrays of the method. This can be used to change internal values as needed.
"""
full_cache(i::DEIntegrator) = error("full_cache: method has not been implemented for the integrator")

"""
    resize_non_user_cache!(integrator::DEIntegrator,k::Int)

Resizes the non-user facing caches to be compatible with a DE of size `k`. This includes resizing Jacobian caches.

!!! note
    In many cases, [`resize!`](@ref) simply resizes [`full_cache`](@ref) variables and then
    calls this function. This finer control is required for some `AbstractArray`
    operations.
"""
resize_non_user_cache!(i::DEIntegrator,ii::Int) = error("resize_non_user_cache!: method has not been implemented for the integrator")

"""
    deleteat_non_user_cache!(integrator::DEIntegrator,idxs)

[`deleteat!`](@ref)s the non-user facing caches at indices `idxs`. This includes resizing Jacobian caches.

!!! note
    In many cases, `deleteat!` simply `deleteat!`s [`full_cache`](@ref) variables and then
    calls this function. This finer control is required for some `AbstractArray`
    operations.
"""
deleteat_non_user_cache!(i::DEIntegrator,idxs) = error("deleteat_non_user_cache!: method has not been implemented for the integrator")

"""
    addat_non_user_cache!(i::DEIntegrator,idxs)

[`addat!`](@ref)s the non-user facing caches at indices `idxs`. This includes resizing Jacobian caches.
!!! note
    In many cases, `addat!` simply `addat!`s [`full_cache`](@ref) variables and then
    calls this function. This finer control is required for some `AbstractArray`
    operations.
"""
addat_non_user_cache!(i::DEIntegrator,idxs) = error("addat_non_user_cache!: method has not been implemented for the integrator")

"""
    terminate!(i::DEIntegrator[, retcode = :Terminated])

Terminates the integrator by emptying `tstops`. This can be used in events and callbacks to immediately
end the solution process.  Optionally, `retcode` may be specified (see: [Return Codes (RetCodes)](@ref retcodes)).
"""
terminate!(i::DEIntegrator) = error("terminate!: method has not been implemented for the integrator")

"""
    get_du(i::DEIntegrator)

Returns the derivative at `t`.
"""
get_du(i::DEIntegrator) = error("get_du: method has not been implemented for the integrator")

"""
    get_du!(out,i::DEIntegrator)

Write the current derivative at `t` into `out`.
"""
get_du!(out,i::DEIntegrator) = error("get_du: method has not been implemented for the integrator")
get_dt(i::DEIntegrator) = error("get_dt: method has not been implemented for the integrator")

"""
    get_proposed_dt(i::DEIntegrator)

Gets the proposed `dt` for the next timestep.
"""
get_proposed_dt(i::DEIntegrator) = error("get_proposed_dt: method has not been implemented for the integrator")

"""
    set_proposed_dt(i::DEIntegrator,dt)
    set_proposed_dt(i::DEIntegrator,i2::DEIntegrator)

Sets the proposed `dt` for the next timestep. If second argument isa `DEIntegrator` then it sets the timestepping of
first argument to match that of second one. Note that due to PI control and step acceleration this is more than matching
the factors in most cases.
"""
set_proposed_dt!(i::DEIntegrator,dt) = error("set_proposed_dt!: method has not been implemented for the integrator")

"""
    savevalues!(integrator::DEIntegrator,
      force_save=false) -> Tuple{Bool, Bool}

Try to save the state and time variables at the current time point, or the
`saveat` point by using interpolation when appropriate. It returns a tuple that
is `(saved, savedexactly)`. If `savevalues!` saved value, then `saved` is true,
and if `savevalues!` saved at the current time point, then `savedexactly` is
true.

The saving priority/order is as follows:
  - `save_on`
    - `saveat`
    - `force_save`
    - `save_everystep`
"""
savevalues!(i::DEIntegrator) = error("savevalues!: method has not been implemented for the integrator")

"""
    u_modified!(i::DEIntegrator,bool)

Sets `bool` which states whether a change to `u` occurred, allowing the solver to handle the discontinuity. By default,
this is assumed to be true if a callback is used. This will result in the re-calculation of the derivative at
`t+dt`, which is not necessary if the algorithm is FSAL and `u` does not experience a discontinuous change at the
end of the interval. Thus if `u` is unmodified in a callback, a single call to the derivative calculation can be
eliminated by `u_modified!(integrator,false)`.
"""
u_modified!(i::DEIntegrator,bool) = error("u_modified!: method has not been implemented for the integrator")

"""
    add_tstop!(i::DEIntegrator,t)

Adds a `tstop` at time `t`.
"""
add_tstop!(i::DEIntegrator,t) = error("add_tstop!: method has not been implemented for the integrator")

"""
    add_saveat!(i::DEIntegrator,t)

Adds a `saveat` time point at `t`.
"""
add_saveat!(i::DEIntegrator,t) = error("add_saveat!: method has not been implemented for the integrator")

set_abstol!(i::DEIntegrator,t) = error("set_abstol!: method has not been implemented for the integrator")
set_reltol!(i::DEIntegrator,t) = error("set_reltol!: method has not been implemented for the integrator")

"""
    reinit!(integrator::DEIntegrator,args...; kwargs...)

The reinit function lets you restart the integration at a new value.

# Arguments

- `u0`: Value of `u` to start at. Default value is `integrator.sol.prob.u0`

# Keyword Arguments
- `t0`: Starting timepoint. Default value is `integrator.sol.prob.tspan[1]`
- `tf`: Ending timepoint. Default value is `integrator.sol.prob.tspan[2]`
- `erase_sol=true`: Whether to start with no other values in the solution, or keep the previous solution.
- `tstops`, `d_discontinuities`, & `saveat`: Cache where these are stored. Default is the original cache.
- `reset_dt`: Set whether to reset the current value of `dt` using the automatic `dt` determination algorithm. Default is
  `(integrator.dtcache == zero(integrator.dt)) && integrator.opts.adaptive`
- `reinit_callbacks`: Set whether to run the callback initializations again (and `initialize_save` is for that). Default is `true`.
- `reinit_cache`: Set whether to re-run the cache initialization function (i.e. resetting FSAL, not allocating vectors)
  which should usually be true for correctness. Default is `true`.

Additionally, once can access [`auto_dt_reset!`](@ref) which will run the auto `dt` initialization algorithm.
"""
reinit!(integrator::DEIntegrator,args...; kwargs...) =
       error("reinit!: method has not been implemented for the integrator")

"""
initialize_dae!(integrator::DEIntegrator,initializealg = integrator.initializealg)

Runs the DAE initialization to find a consistent state vector. The optional
argument `initializealg` can be used to specify a different initialization
algorithm to use.
"""
initialize_dae!(integrator::DEIntegrator) =
       error("initialize_dae!: method has not been implemented for the integrator")

"""
    auto_dt_reset!(integrator::DEIntegrator)

Run the auto `dt` initialization algorithm.
"""
auto_dt_reset!(integrator::DEIntegrator) = error("auto_dt_reset!: method has not been implemented for the integrator")

"""
    change_t_via_interpolation!(integrator::DEIntegrator,t,modify_save_endpoint=Val{false})

Modifies the current `t` and changes all of the corresponding values using the local interpolation. If the current solution
has already been saved, one can provide the optional value `modify_save_endpoint` to also modify the endpoint of `sol` in the
same manner.
"""
change_t_via_interpolation!(i::DEIntegrator,args...) = error("change_t_via_interpolation!: method has not been implemented for the integrator")

addsteps!(i::DEIntegrator,args...) = nothing

"""
    reeval_internals_due_to_modification!(integrator::DDEIntegrator)

Recalculate interpolation data and update ODE integrator after changes by callbacks.
"""
reeval_internals_due_to_modification!(integrator::DEIntegrator) = nothing

"""
    set_t!(integrator::DEIntegrator, t)

Set current time point of the `integrator` to `t`.
"""
set_t!(integrator::DEIntegrator, t) =
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

function addat!(a::AbstractArray,idxs,val=nothing)
  if val === nothing
    resize!(a,length(a)+length(idxs))
  else
    error("real addat! on arrays isn't supported yet")
    #=
    flip_range = last(idxs):-1:idxs.start
    @show idxs,flip_range
    splice!(a,flip_range,val)
    =#
  end
end

### Integrator traits

has_reinit(i::DEIntegrator) = false

### Display

Base.summary(I::DEIntegrator) = string(
                  TYPE_COLOR, nameof(typeof(I)),
                  NO_COLOR, " with uType ",
                  TYPE_COLOR, typeof(I.u),
                  NO_COLOR, " and tType ",
                  TYPE_COLOR, typeof(I.t), NO_COLOR)

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
TreeViews.hastreeview(x::DiffEqBase.DEIntegrator) = true
function TreeViews.treelabel(io::IO,x::DiffEqBase.DEIntegrator,
                             mime::MIME"text/plain" = MIME"text/plain"())
  show(io,mime,Base.Text(Base.summary(x)))
end

### Error check (retcode)

last_step_failed(integrator::DEIntegrator) = false

"""
    check_error(integrator)

Check state of `integrator` and return one of the
[Return Codes](http://docs.juliadiffeq.org/dev/basics/solution.html#Return-Codes-(RetCodes)-1)
"""
function check_error(integrator::DEIntegrator)
  # This implementation is intended to be used for ODEIntegrator and
  # SDEIntegrator.
  if isnan(integrator.dt)
    if integrator.opts.verbose
      @warn("NaN dt detected. Likely a NaN value in the state, parameters, or derivative value caused this outcome.")
    end
    return :DtNaN
  end
  if integrator.iter > integrator.opts.maxiters
    if integrator.opts.verbose
      @warn("Interrupted. Larger maxiters is needed.")
    end
    return :MaxIters
  end

  # The last part:
  # If you are close to the end, don't exit: let the user hit the end!
  # However, if we try that and the step fails, exit instead of infinite loop
  if !integrator.opts.force_dtmin && integrator.opts.adaptive &&
     abs(integrator.dt) <= abs(integrator.opts.dtmin) &&
     (((hasproperty(integrator,:opts) && hasproperty(integrator.opts,:tstops)) ?
     integrator.t + integrator.dt < integrator.tdir*first(integrator.opts.tstops) :
     true) || (hasproperty(integrator,:accept_step) && !integrator.accept_step))
    if integrator.opts.verbose
      @warn("dt <= dtmin. Aborting. There is either an error in your model specification or the true solution is unstable.")
    end
    return :DtLessThanMin
  end
  if integrator.opts.unstable_check(integrator.dt,integrator.u,integrator.p,integrator.t)
    if integrator.opts.verbose
      @warn("Instability detected. Aborting")
    end
    return :Unstable
  end
  if last_step_failed(integrator)
    if integrator.opts.verbose
      @warn("Newton steps could not converge and algorithm is not adaptive. Use a lower dt.")
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
function done(integrator::DEIntegrator)
  if ! (integrator.sol.retcode in (:Default, :Success))
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
function Base.iterate(integrator::DEIntegrator,state=0)
  done(integrator) && return nothing
  state += 1
  step!(integrator) # Iter updated in the step! header
  # Next is callbacks -> iterator  -> top
  return integrator,state
end

Base.eltype(::Type{T}) where {T<:DEIntegrator} = T
Base.IteratorSize(::Type{<:DEIntegrator}) = Base.SizeUnknown()

### Other Iterators

struct IntegratorTuples{I}
 integrator::I
end

function Base.iterate(tup::IntegratorTuples, state=0)
  done(tup.integrator) && return nothing
  step!(tup.integrator) # Iter updated in the step! header
  state += 1
  # Next is callbacks -> iterator  -> top
  return (tup.integrator.u,tup.integrator.t),state
end

Base.eltype(::Type{IntegratorTuples{I}}) where {U, T, I<:DEIntegrator{<:Any, <:Any, U, T}} = Tuple{U, T} 
Base.IteratorSize(::Type{<:IntegratorTuples}) = Base.SizeUnknown()

RecursiveArrayTools.tuples(integrator::DEIntegrator) = IntegratorTuples(integrator)

"""
$(TYPEDEF)
"""
struct IntegratorIntervals{I}
 integrator::I
end

function Base.iterate(tup::IntegratorIntervals,state=0)
  done(tup.integrator) && return nothing
  state += 1
  step!(tup.integrator) # Iter updated in the step! header
  # Next is callbacks -> iterator  -> top
  return (tup.integrator.uprev,tup.integrator.tprev,tup.integrator.u,tup.integrator.t),state
end

Base.eltype(::Type{IntegratorIntervals{I}}) where {U, T, I<:DEIntegrator{<:Any, <:Any, U, T}} = Tuple{U, T, U, T} 
Base.IteratorSize(::Type{<:IntegratorIntervals}) = Base.SizeUnknown()

intervals(integrator::DEIntegrator) = IntegratorIntervals(integrator)

struct TimeChoiceIterator{T,T2}
  integrator::T
  ts::T2
end

function Base.iterate(iter::TimeChoiceIterator,state=1)
  state > length(iter.ts) && return nothing
  t = iter.ts[state]
  integrator = iter.integrator
  if isinplace(integrator.sol.prob)
    tmp = first(get_tmp_cache(integrator))
    if t == integrator.t
      tmp .= integrator.u
    elseif t < integrator.t
      integrator(tmp,t)
    else
      step!(integrator,t-integrator.t)
      integrator(tmp,t)
    end
    return (tmp,t),state+1
  else
    if t == integrator.t
      tmp = integrator.u
    elseif t < integrator.t
      tmp = integrator(t)
    else
      step!(integrator,t-integrator.t)
      tmp = integrator(t)
    end
    return (tmp,t),state+1
  end
end

Base.length(iter::TimeChoiceIterator) = length(iter.ts)

@recipe function f(integrator::DEIntegrator;
                    denseplot=(integrator.opts.calck || typeof(integrator) <: AbstractSDEIntegrator)  && integrator.iter>0,
                    plotdensity =10,
                    plot_analytic=false,vars=nothing)

  int_vars = interpret_vars(vars,integrator.sol)

  if denseplot
    # Generate the points from the plot from dense function
    plott = collect(range(integrator.tprev;step=integrator.t,length=plotdensity))
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
  (plot_vecs...,)
end

function step!(integ::DEIntegrator, dt, stop_at_tdt = false)
    (dt * integ.tdir) < 0 * oneunit(dt) && error("Cannot step backward.")
    t = integ.t
    next_t = t+dt
    stop_at_tdt && add_tstop!(integ,next_t)
    while integ.t*integ.tdir < next_t*integ.tdir
        step!(integ)
        integ.sol.retcode in (:Default, :Success) || break
    end
end

has_destats(i::DEIntegrator) = false

"""
    is_integrator_adaptive(i::DEIntegrator)

Checks if the integrator is adaptive
"""
isadaptive(integrator::DEIntegrator) =
        isdefined(integrator.opts, :adaptive) ? integrator.opts.adaptive : false
