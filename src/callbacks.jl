# Necessary to have initialize set u_modified to false if all don't do anything
# otherwise unnecessary save
INITIALIZE_DEFAULT(cb,u,t,integrator) = u_modified!(integrator, false)
FINALIZE_DEFAULT(cb,u,t,integrator) = nothing

@enum RootfindOpt::Int8 begin
  NoRootFind    = 0
  LeftRootFind  = 1
  RightRootFind = 2
end

function Base.convert(::Type{RootfindOpt}, b::Bool)
  return b ? LeftRootFind : NoRootFind
end

"""
```julia
ContinuousCallback(condition,affect!,affect_neg!;
                   initialize = INITIALIZE_DEFAULT,
                   finalize = FINALIZE_DEFAULT,
                   idxs = nothing,
                   rootfind=LeftRootFind,
                   save_positions=(true,true),
                   interp_points=10,
                   abstol=10eps(),reltol=0)
```

```julia
function ContinuousCallback(condition,affect!;
                   initialize = INITIALIZE_DEFAULT,
                   finalize = FINALIZE_DEFAULT,
                   idxs = nothing,
                   rootfind=LeftRootFind,
                   save_positions=(true,true),
                   affect_neg! = affect!,
                   interp_points=10,
                   abstol=10eps(),reltol=0)
```

Contains a single callback whose `condition` is a continuous function. The callback is triggered when this function evaluates to 0.

# Arguments
- `condition`: This is a function `condition(u,t,integrator)` for declaring when
  the callback should be used. A callback is initiated if the condition hits
  `0` within the time interval. See the [Integrator Interface](@ref integrator) documentation for information about `integrator`.
- `affect!`: This is the function `affect!(integrator)` where one is allowed to
  modify the current state of the integrator. If you do not pass an `affect_neg!`
  function, it is called when `condition` is found to be `0` (at a root) and
  the cross is either an upcrossing (from negative to positive) or a downcrossing
  (from positive to negative). You need to explicitly pass `nothing` as the
  `affect_neg!` argument if it should only be called at upcrossings, e.g.
  `ContinuousCallback(condition, affect!, nothing)`. For more information on what can
  be done, see the [Integrator Interface](@ref integrator) manual page. Modifications to
  `u` are safe in this function.
- `affect_neg!=affect!`: This is the function `affect_neg!(integrator)` where one is allowed to
  modify the current state of the integrator. This is called when `condition` is
  found to be `0` (at a root) and the cross is an downcrossing (from positive to
  negative). For more information on what can
  be done, see the [Integrator Interface](@ref integrator) manual page. Modifications to
  `u` are safe in this function.
- `rootfind=LeftRootFind`: This is a flag to specify the type of rootfinding to do for finding
  event location. If this is set to `LeftRootfind`, the solution will be backtracked to the point where
  `condition==0` and if the solution isn't exact, the left limit of root is used. If set to
  `RightRootFind`, the solution would be set to the right limit of the root. Otherwise the systems and
  the `affect!` will occur at `t+dt`.
- `interp_points=10`: The number of interpolated points to check the condition. The
  condition is found by checking whether any interpolation point / endpoint has
  a different sign. If `interp_points=0`, then conditions will only be noticed if
  the sign of `condition` is different at `t` than at `t+dt`. This behavior is not
  robust when the solution is oscillatory, and thus it's recommended that one use
  some interpolation points (they're cheap to compute!).
  `0` within the time interval.
- `save_positions=(true,true)`: Boolean tuple for whether to save before and after the `affect!`.
  This saving will occur just before and after the event, only at event times, and
  does not depend on options like `saveat`, `save_everystep`, etc. (i.e. if
  `saveat=[1.0,2.0,3.0]`, this can still add a save point at `2.1` if true).
  For discontinuous changes like a modification to `u` to be
  handled correctly (without error), one should set `save_positions=(true,true)`.
- `idxs=nothing`: The components which will be interpolated into the condition. Defaults
  to `nothing` which means `u` will be all components.
- `initialize`: This is a function `(c,u,t,integrator)` which can be used to initialize
  the state of the callback `c`. It should modify the argument `c` and the return is
  ignored.
- `finalize`: This is a function `(c,u,t,integrator)` which can be used to finalize
  the state of the callback `c`. It can modify the argument `c` and the return is ignored.
- `abstol=1e-14` & `reltol=0`: These are used to specify a tolerance from zero for the rootfinder:
  if the starting condition is less than the tolerance from zero, then no root will be detected.
  This is to stop repeat events happening just after a previously rootfound event.
"""
struct ContinuousCallback{F1,F2,F3,F4,F5,T,T2,I,R} <: AbstractContinuousCallback
  condition::F1
  affect!::F2
  affect_neg!::F3
  initialize::F4
  finalize::F5
  idxs::I
  rootfind::RootfindOpt
  interp_points::Int
  save_positions::BitArray{1}
  dtrelax::R
  abstol::T
  reltol::T2
  ContinuousCallback(condition::F1,affect!::F2,affect_neg!::F3,
                     initialize::F4,finalize::F5,idxs::I,rootfind,
                     interp_points,save_positions,dtrelax::R,abstol::T,reltol::T2) where {F1,F2,F3,F4,F5,T,T2,I,R} =
                       new{F1,F2,F3,F4,F5,T,T2,I,R}(condition,
                                               affect!,affect_neg!,
                                               initialize,finalize,idxs,rootfind,interp_points,
                                               BitArray(collect(save_positions)),
                                               dtrelax,abstol,reltol)
end

ContinuousCallback(condition,affect!,affect_neg!;
                   initialize = INITIALIZE_DEFAULT,
                   finalize = FINALIZE_DEFAULT,
                   idxs = nothing,
                   rootfind=LeftRootFind,
                   save_positions=(true,true),
                   interp_points=10,
                   dtrelax=1,
                   abstol=10eps(),reltol=0) = ContinuousCallback(
                              condition,affect!,affect_neg!,initialize,finalize,
                              idxs,
                              rootfind,interp_points,
                              save_positions,
                              dtrelax,abstol,reltol)

function ContinuousCallback(condition,affect!;
                   initialize = INITIALIZE_DEFAULT,
                   finalize = FINALIZE_DEFAULT,
                   idxs = nothing,
                   rootfind=LeftRootFind,
                   save_positions=(true,true),
                   affect_neg! = affect!,
                   interp_points=10,
                   dtrelax=1,
                   abstol=10eps(),reltol=0)

 ContinuousCallback(
            condition,affect!,affect_neg!,initialize,finalize,idxs,
            rootfind,interp_points,
            collect(save_positions),
            dtrelax,abstol,reltol)

end

"""
```julia
VectorContinuousCallback(condition,affect!,affect_neg!,len;
                         initialize = INITIALIZE_DEFAULT,
                         finalize = FINALIZE_DEFAULT,
                         idxs = nothing,
                         rootfind=LeftRootFind,
                         save_positions=(true,true),
                         interp_points=10,
                         abstol=10eps(),reltol=0)
```

```julia
VectorContinuousCallback(condition,affect!,len;
                   initialize = INITIALIZE_DEFAULT,
                   finalize = FINALIZE_DEFAULT,
                   idxs = nothing,
                   rootfind=LeftRootFind,
                   save_positions=(true,true),
                   affect_neg! = affect!,
                   interp_points=10,
                   abstol=10eps(),reltol=0)
```

This is also a subtype of `AbstractContinuousCallback`. `CallbackSet` is not feasible when you have a large number of callbacks,
as it doesn't scale well. For this reason, we have `VectorContinuousCallback` - it allows you to have a single callback for
multiple events.

# Arguments

- `condition`: This is a function `condition(out, u, t, integrator)` which should save the condition value in the array `out`
   at the right index. Maximum index of `out` should be specified in the `len` property of callback. So this way you can have
   a chain of `len` events, which would cause the `i`th event to trigger when `out[i] = 0`.
- `affect!`: This is a function `affect!(integrator, event_index)` which lets you modify `integrator` and it tells you about
   which event occured using `event_idx` i.e. gives you index `i` for which `out[i]` came out to be zero.
- `len`: Number of callbacks chained. This is compulsory to be specified.

Rest of the arguments have the same meaning as in [`ContinuousCallback`](@ref).
"""
struct VectorContinuousCallback{F1,F2,F3,F4,F5,T,T2,I,R} <: AbstractContinuousCallback
  condition::F1
  affect!::F2
  affect_neg!::F3
  len::Int
  initialize::F4
  finalize::F5
  idxs::I
  rootfind::RootfindOpt
  interp_points::Int
  save_positions::BitArray{1}
  dtrelax::R
  abstol::T
  reltol::T2
  VectorContinuousCallback(condition::F1,affect!::F2,affect_neg!::F3,len::Int,
                           initialize::F4,finalize::F5,idxs::I,rootfind,
                           interp_points,save_positions,dtrelax::R,
                           abstol::T,reltol::T2) where {F1,F2,F3,F4,F5,T,T2,I,R} =
                       new{F1,F2,F3,F4,F5,T,T2,I,R}(condition,
                                               affect!,affect_neg!,len,
                                               initialize,finalize,idxs,rootfind,interp_points,
                                               BitArray(collect(save_positions)),
                                               dtrelax,abstol,reltol)
end

VectorContinuousCallback(condition,affect!,affect_neg!,len;
                         initialize = INITIALIZE_DEFAULT,
                         finalize = FINALIZE_DEFAULT,
                         idxs = nothing,
                         rootfind=LeftRootFind,
                         save_positions=(true,true),
                         interp_points=10,
                         dtrelax=1,
                         abstol=10eps(),reltol=0) = VectorContinuousCallback(
                              condition,affect!,affect_neg!,len,
                              initialize,finalize,
                              idxs,
                              rootfind,interp_points,
                              save_positions,dtrelax,
                              abstol,reltol)

function VectorContinuousCallback(condition,affect!,len;
                   initialize = INITIALIZE_DEFAULT,
                   finalize = FINALIZE_DEFAULT,
                   idxs = nothing,
                   rootfind=LeftRootFind,
                   save_positions=(true,true),
                   affect_neg! = affect!,
                   interp_points=10,
                   dtrelax=1,
                   abstol=10eps(),reltol=0)

 VectorContinuousCallback(
            condition,affect!,affect_neg!,len,initialize,finalize,idxs,
            rootfind,interp_points,
            collect(save_positions),
            dtrelax,abstol,reltol)

end

"""
```julia
DiscreteCallback(condition,affect!;
                 initialize = INITIALIZE_DEFAULT,
                 finalize = FINALIZE_DEFAULT,
                 save_positions=(true,true))
```

# Arguments

- `condition`: This is a function `condition(u,t,integrator)` for declaring when
  the callback should be used. A callback is initiated if the condition evaluates
  to `true`. See the [Integrator Interface](@ref integrator) documentation for information about `integrator`.
    - `affect!`: This is the function `affect!(integrator)` where one is allowed to
  modify the current state of the integrator. For more information on what can
  be done, see the [Integrator Interface](@ref integrator) manual page.
- `save_positions`: Boolean tuple for whether to save before and after the `affect!`.
  This saving will occur just before and after the event, only at event times, and
  does not depend on options like `saveat`, `save_everystep`, etc. (i.e. if
  `saveat=[1.0,2.0,3.0]`, this can still add a save point at `2.1` if true).
  For discontinuous changes like a modification to `u` to be
  handled correctly (without error), one should set `save_positions=(true,true)`.
- `initialize`: This is a function `(c,u,t,integrator)` which can be used to initialize
  the state of the callback `c`. It should modify the argument `c` and the return is
  ignored.
- `finalize`: This is a function `(c,u,t,integrator)` which can be used to finalize
  the state of the callback `c`. It should can the argument `c` and the return is
  ignored.
"""
struct DiscreteCallback{F1,F2,F3,F4} <: AbstractDiscreteCallback
  condition::F1
  affect!::F2
  initialize::F3
  finalize::F4
  save_positions::BitArray{1}
  DiscreteCallback(condition::F1,affect!::F2,
                   initialize::F3,finalize::F4,save_positions) where {F1,F2,F3,F4} = new{F1,F2,F3,F4}(condition,
                                                                                   affect!,initialize,finalize,
                                                                                   BitArray(collect(save_positions)))
end
DiscreteCallback(condition,affect!;
                 initialize = INITIALIZE_DEFAULT, finalize = FINALIZE_DEFAULT,
                 save_positions=(true,true)) = DiscreteCallback(condition,affect!,initialize,finalize,save_positions)

"""
$(TYPEDEF)

Multiple callbacks can be chained together to form a `CallbackSet`. A `CallbackSet`
is constructed by passing the constructor `ContinuousCallback`, `DiscreteCallback`,
`VectorContinuousCallback` or other `CallbackSet` instances:

    CallbackSet(cb1,cb2,cb3)

You can pass as many callbacks as you like. When the solvers encounter multiple
callbacks, the following rules apply:

* `ContinuousCallback`s and `VectorContinuousCallback`s are applied before `DiscreteCallback`s. (This is because
  they often implement event-finding that will backtrack the timestep to smaller
  than `dt`).
* For `ContinuousCallback`s and `VectorContinuousCallback`s, the event times are found by rootfinding and only
  the first `ContinuousCallback` or `VectorContinuousCallback` affect is applied.
* The `DiscreteCallback`s are then applied in order. Note that the ordering only
  matters for the conditions: if a previous callback modifies `u` in such a way
  that the next callback no longer evaluates condition to `true`, its `affect`
  will not be applied.

"""
struct CallbackSet{T1<:Tuple,T2<:Tuple} <: DECallback
  continuous_callbacks::T1
  discrete_callbacks::T2
end

CallbackSet(callback::AbstractDiscreteCallback) = CallbackSet((),(callback,))
CallbackSet(callback::AbstractContinuousCallback) = CallbackSet((callback,),())
CallbackSet() = CallbackSet((),())
CallbackSet(cb::Nothing) = CallbackSet()

# For Varargs, use recursion to make it type-stable
CallbackSet(callbacks::Union{DECallback,Nothing}...) = CallbackSet(split_callbacks((), (), callbacks...)...)

"""
    split_callbacks(cs, ds, args...)

Split comma seperated callbacks into sets of continous and discrete callbacks.
"""
@inline split_callbacks(cs, ds) = cs, ds
@inline split_callbacks(cs, ds, c::Nothing, args...) = split_callbacks(cs, ds, args...)
@inline split_callbacks(cs, ds, c::AbstractContinuousCallback, args...) = split_callbacks((cs..., c), ds, args...)
@inline split_callbacks(cs, ds, d::AbstractDiscreteCallback, args...) = split_callbacks(cs, (ds..., d), args...)
@inline function split_callbacks(cs, ds, d::CallbackSet, args...)
  split_callbacks((cs...,d.continuous_callbacks...), (ds..., d.discrete_callbacks...), args...)
end

"""
    initialize!(cb::CallbackSet,u,t,integrator::DEIntegrator)

Recursively apply `initialize!` and return whether any modified u
"""
function initialize!(cb::CallbackSet,u,t,integrator::DEIntegrator)
  initialize!(u,t,integrator,false,cb.continuous_callbacks...,cb.discrete_callbacks...)
end
initialize!(cb::CallbackSet{Tuple{},Tuple{}},u,t,integrator::DEIntegrator) = false
function initialize!(u,t,integrator::DEIntegrator,any_modified::Bool,
                     c::DECallback,cs::DECallback...)
  c.initialize(c,u,t,integrator)
  initialize!(u,t,integrator,any_modified || integrator.u_modified,cs...)
end
function initialize!(u,t,integrator::DEIntegrator,any_modified::Bool,
                     c::DECallback)
  c.initialize(c,u,t,integrator)
  any_modified || integrator.u_modified
end


"""
    finalize!(cb::CallbackSet,u,t,integrator::DEIntegrator)

Recursively apply `finalize!` and return whether any modified u
"""
function finalize!(cb::CallbackSet,u,t,integrator::DEIntegrator)
  finalize!(u,t,integrator,false,cb.continuous_callbacks...,cb.discrete_callbacks...)
end
finalize!(cb::CallbackSet{Tuple{},Tuple{}},u,t,integrator::DEIntegrator) = false
function finalize!(u,t,integrator::DEIntegrator,any_modified::Bool,
                     c::DECallback,cs::DECallback...)
  c.finalize(c,u,t,integrator)
  finalize!(u,t,integrator,any_modified || integrator.u_modified,cs...)
end
function finalize!(u,t,integrator::DEIntegrator,any_modified::Bool,
                     c::DECallback)
  c.finalize(c,u,t,integrator)
  any_modified || integrator.u_modified
end


# Helpers
Base.isempty(cb::CallbackSet) = isempty(cb.continuous_callbacks) && isempty(cb.discrete_callbacks)
Base.isempty(cb::AbstractContinuousCallback) = false
Base.isempty(cb::AbstractDiscreteCallback) = false

has_continuous_callback(cb::DiscreteCallback) = false
has_continuous_callback(cb::ContinuousCallback) = true
has_continuous_callback(cb::VectorContinuousCallback) = true
has_continuous_callback(cb::CallbackSet) = !isempty(cb.continuous_callbacks)
has_continuous_callback(cb::Nothing) = false

#======================================================#
# Callback handling
#======================================================#

function get_tmp(integrator::DEIntegrator, callback)
  _tmp = get_tmp_cache(integrator)
  _tmp === nothing && return nothing
  _cache = first(_tmp)
  if callback.idxs === nothing
    tmp = _cache
  elseif !(callback.idxs isa Number)
    tmp = @view _cache[callback.idxs]
  else
    tmp = nothing
  end
  return tmp
end

function get_condition(integrator::DEIntegrator, callback, abst)
  tmp = get_tmp(integrator, callback)
  ismutable = !(tmp === nothing)
  if abst == integrator.t
    if callback.idxs === nothing
      tmp = integrator.u
    elseif callback.idxs isa Number
      tmp = integrator.u[callback.idxs]
    else
      tmp = @view integrator.u[callback.idxs]
    end
  elseif abst == integrator.tprev
    if callback.idxs === nothing
      tmp = integrator.uprev
    elseif callback.idxs isa Number
      tmp = integrator.uprev[callback.idxs]
    else
      tmp = @view integrator.uprev[callback.idxs]
    end
  else
    if ismutable
      if callback.idxs === nothing
        integrator(tmp,abst,Val{0})
      else
        integrator(tmp,abst,Val{0},idxs=callback.idxs)
      end
    else
      if callback.idxs === nothing
        tmp = integrator(abst,Val{0})
      else
        tmp = integrator(abst,Val{0},idxs=callback.idxs)
      end
    end
    # ismutable && !(callback.idxs isa Number) ? integrator(tmp,abst,Val{0},idxs=callback.idxs) :
    #                                                 tmp = integrator(abst,Val{0},idxs=callback.idxs)
  end
  integrator.sol.destats.ncondition += 1
  if callback isa VectorContinuousCallback
    callback.condition(@view(integrator.callback_cache.tmp_condition[1:callback.len]),tmp,abst,integrator)
    return @view(integrator.callback_cache.tmp_condition[1:callback.len])
  else
    return callback.condition(tmp,abst,integrator)
  end
end

# Use Recursion to find the first callback for type-stability

# Base Case: Only one callback
function find_first_continuous_callback(integrator, callback::AbstractContinuousCallback)
  (find_callback_time(integrator,callback,1)...,1,1)
end

# Starting Case: Compute on the first callback
function find_first_continuous_callback(integrator, callback::AbstractContinuousCallback, args...)
  find_first_continuous_callback(integrator,find_callback_time(integrator,callback,1)...,1,1,args...)
end

function find_first_continuous_callback(integrator,tmin::Number,upcrossing::Number,
                                        event_occured::Bool,event_idx::Int,idx::Int,counter::Int,
                                        callback2)
  counter += 1 # counter is idx for callback2.
  tmin2,upcrossing2,event_occurred2,event_idx2 = find_callback_time(integrator,callback2,counter)

  if event_occurred2 && (tmin2 < tmin || !event_occured)
    return tmin2,upcrossing2,true,event_idx2,counter,counter
  else
    return tmin,upcrossing,event_occured,event_idx,idx,counter
  end
end

function find_first_continuous_callback(integrator,tmin::Number,upcrossing::Number,event_occured::Bool,event_idx::Int,idx::Int,counter::Int,callback2,args...)
  find_first_continuous_callback(integrator,find_first_continuous_callback(integrator,tmin,upcrossing,event_occured,event_idx,idx,counter,callback2)...,args...)
end

@inline function determine_event_occurance(integrator,callback::VectorContinuousCallback,counter)
  event_occurred = false
  if callback.interp_points!=0
    addsteps!(integrator)
  end
  ts = range(integrator.tprev, stop=integrator.t, length=callback.interp_points)
  interp_index = 0
  # Check if the event occured
  previous_condition = @views(integrator.callback_cache.previous_condition[1:callback.len])

  if callback.idxs === nothing
    callback.condition(previous_condition,integrator.uprev,integrator.tprev,integrator)
  else
    callback.condition(previous_condition,integrator.uprev[callback.idxs],integrator.tprev,integrator)
  end
  integrator.sol.destats.ncondition += 1

  ivec = integrator.vector_event_last_time
  prev_sign = @view(integrator.callback_cache.prev_sign[1:callback.len])
  next_sign = @view(integrator.callback_cache.next_sign[1:callback.len])

  if integrator.event_last_time == counter && minimum(ODE_DEFAULT_NORM(ArrayInterface.allowed_getindex(previous_condition,ivec),integrator.t)) <= 100ODE_DEFAULT_NORM(integrator.last_event_error,integrator.t)

    # If there was a previous event, utilize the derivative at the start to
    # chose the previous sign. If the derivative is positive at tprev, then
    # we treat `prev_sign` as negetive, and if the derivative is negative then we
    # treat `prev_sign` as positive, regardless of the postiivity/negativity
    # of the true value due to it being =0 sans floating point issues.

    # Only due this if the discontinuity did not move it far away from an event
    # Since near even we use direction instead of location to reset

    if callback.interp_points==0
      addsteps!(integrator)
    end

    # Evaluate condition slightly in future
    if integrator.t == 0
      abst = integrator.tprev+integrator.tdir*abs(integrator.dt/10000)
    else
      abst = integrator.tprev+integrator.tdir*100*eps(integrator.t)
    end
    tmp_condition = get_condition(integrator, callback, abst)

    # Sometimes users may "switch off" the condition after crossing
    # This is necessary to ensure proper non-detection of a root
    # == is for exact floating point equality!
    @. prev_sign = sign(previous_condition)
    prev_sign[ivec] = tmp_condition[ivec] > previous_condition[ivec] ? 1.0 :
                  (tmp_condition[ivec] == previous_condition[ivec] ?
                  (prev_sign[ivec] = sign(previous_condition[ivec])) : -1.0)
  else
      @. prev_sign = sign(previous_condition)
  end

  prev_sign_index = 1
  abst = integrator.t
  next_condition = get_condition(integrator, callback, abst)
  @. next_sign = sign(next_condition)

  event_idx = findall_events(callback.affect!,callback.affect_neg!,prev_sign,next_sign)
  if length(event_idx) != 0
    event_occurred = true
    interp_index = callback.interp_points
  end
  if callback.interp_points!=0 && !isdiscrete(integrator.alg) && length(prev_sign) != length(event_idx) # Use the interpolants for safety checking
    for i in 2:length(ts)
      abst = ts[i]
      new_sign = get_condition(integrator, callback, abst)
      _event_idx = findall_events(callback.affect!,callback.affect_neg!,prev_sign,new_sign)
      if length(_event_idx) != 0
        event_occurred = true
        event_idx = _event_idx
        interp_index = i
        break
      else
        prev_sign_index = i
      end
    end
  end

  event_idx_out = convert(Array,event_idx) # No-op on arrays
  event_occurred,interp_index,ts,prev_sign,prev_sign_index,event_idx_out
end

@inline function determine_event_occurance(integrator,callback::ContinuousCallback,counter)
  event_occurred = false
  if callback.interp_points!=0
    addsteps!(integrator)
  end
  ts = range(integrator.tprev, stop=integrator.t, length=callback.interp_points)
  interp_index = 0
  # Check if the event occured
  if callback.idxs === nothing
    previous_condition = callback.condition(integrator.uprev,integrator.tprev,integrator)
  else
    @views previous_condition = callback.condition(integrator.uprev[callback.idxs],integrator.tprev,integrator)
  end
  integrator.sol.destats.ncondition += 1

  prev_sign = 0.0
  next_sign = 0.0

  if integrator.event_last_time == counter && minimum(ODE_DEFAULT_NORM(previous_condition,integrator.t)) <= 100ODE_DEFAULT_NORM(integrator.last_event_error,integrator.t)

    # If there was a previous event, utilize the derivative at the start to
    # chose the previous sign. If the derivative is positive at tprev, then
    # we treat `prev_sign` as negetive, and if the derivative is negative then we
    # treat `prev_sign` as positive, regardless of the postiivity/negativity
    # of the true value due to it being =0 sans floating point issues.

    # Only due this if the discontinuity did not move it far away from an event
    # Since near even we use direction instead of location to reset

    if callback.interp_points==0
      addsteps!(integrator)
    end

    # Evaluate condition slightly in future
    if integrator.t == 0
      abst = integrator.tprev+integrator.tdir*abs(integrator.dt/10000)
    else
      abst = integrator.tprev+integrator.tdir*100*eps(integrator.t)
    end
    tmp_condition = get_condition(integrator, callback, abst)

    # Sometimes users may "switch off" the condition after crossing
    # This is necessary to ensure proper non-detection of a root
    # == is for exact floating point equality!
    prev_sign =    tmp_condition > previous_condition ? 1.0 :
                  (tmp_condition == previous_condition ?
                  (prev_sign = sign(previous_condition)) : -1.0)
  else
    prev_sign = sign(previous_condition)
  end

  prev_sign_index = 1
  abst = integrator.t
  next_condition = get_condition(integrator, callback, abst)
  next_sign = sign(next_condition)

  if ((prev_sign < 0 && callback.affect! !== nothing) || (prev_sign > 0 && callback.affect_neg! !== nothing)) && prev_sign*next_sign<=0
    event_occurred = true
    interp_index = callback.interp_points
  elseif callback.interp_points!=0 && !isdiscrete(integrator.alg) # Use the interpolants for safety checking
    for i in 2:length(ts)
      abst = ts[i]
      new_sign = get_condition(integrator, callback, abst)
      if ((prev_sign < 0 && callback.affect! !== nothing) || (prev_sign > 0 && callback.affect_neg! !== nothing)) && prev_sign*new_sign<0
        event_occurred = true
        interp_index = i
        break
      else
        prev_sign_index = i
      end
    end
  end
  event_idx = 1

  event_occurred,interp_index,ts,prev_sign,prev_sign_index,event_idx
end

# rough implementation, needs multiple type handling
# always ensures that if r = bisection(f, (x0, x1))
# then either f(nextfloat(r)) == 0 or f(nextfloat(r)) * f(r) < 0
function bisection(f, tup, t_forward::Bool, rootfind::RootfindOpt, abstol, reltol; maxiters=1000)
  if rootfind == LeftRootFind
    NonlinearSolve.solve(NonlinearSolve.NonlinearProblem{false}(f, tup), NonlinearSolve.Falsi(), abstol=abstol, reltol=reltol).left
  else
    NonlinearSolve.solve(NonlinearSolve.NonlinearProblem{false}(f, tup), NonlinearSolve.Falsi(), abstol=abstol, reltol=reltol).right
  end
end

## Different definition for GPUs
function findall_events(affect!,affect_neg!,prev_sign,next_sign)
  findall(x-> ((prev_sign[x] < 0 && affect! !== nothing) || (prev_sign[x] > 0 && affect_neg! !== nothing)) && prev_sign[x]*next_sign[x]<=0, keys(prev_sign))
end

function find_callback_time(integrator,callback::ContinuousCallback,counter)
  event_occurred,interp_index,ts,prev_sign,prev_sign_index,event_idx = determine_event_occurance(integrator,callback,counter)
  if event_occurred
    if callback.condition === nothing
      new_t = zero(typeof(integrator.t))
    else
      if callback.interp_points!=0
        top_t = ts[interp_index] # Top at the smallest
        bottom_t = ts[prev_sign_index]
      else
        top_t = integrator.t
        bottom_t = integrator.tprev
      end
      if callback.rootfind != NoRootFind && !isdiscrete(integrator.alg)
        zero_func(abst, p=nothing) = get_condition(integrator, callback, abst)
        if zero_func(top_t) == 0
          Θ = top_t
        else
          if integrator.event_last_time == counter &&
            abs(zero_func(bottom_t)) <= 100abs(integrator.last_event_error) &&
            prev_sign_index == 1

            # Determined that there is an event by derivative
            # But floating point error may make the end point negative

            sign_top = sign(zero_func(top_t))
            diff_t = integrator.tdir*2eps(bottom_t)
            bottom_t += diff_t
            iter = 1
            # This check should match the same check in bisection
            while sign(zero_func(bottom_t)) * sign_top >= zero(sign_top) && iter < 12
              diff_t *= 5
              bottom_t = integrator.tprev + diff_t
              iter += 1
            end
            iter == 12 && error("Double callback crossing floating pointer reducer errored. Report this issue.")
          end
          Θ = bisection(zero_func, (bottom_t, top_t), isone(integrator.tdir), callback.rootfind, callback.abstol, callback.reltol)
          integrator.last_event_error = ODE_DEFAULT_NORM(zero_func(Θ), Θ)
        end
        #Θ = prevfloat(...)
        # prevfloat guerentees that the new time is either 1 floating point
        # numbers just before the event or directly at zero, but not after.
        # If there's a barrier which is never supposed to be crossed,
        # then this will ensure that
        # The item never leaves the domain. Otherwise Roots.jl can return
        # a float which is slightly after, making it out of the domain, causing
        # havoc.
        new_t = Θ - integrator.tprev
      elseif interp_index != callback.interp_points && !isdiscrete(integrator.alg)
        new_t = ts[interp_index] - integrator.tprev
      else
        # If no solve and no interpolants, just use endpoint
        new_t = integrator.dt
      end
    end
  else
    new_t = zero(typeof(integrator.t))
  end

  new_t,prev_sign,event_occurred,event_idx
end

function find_callback_time(integrator,callback::VectorContinuousCallback,counter)
  event_occurred,interp_index,ts,prev_sign,prev_sign_index,event_idx = determine_event_occurance(integrator,callback,counter)
  if event_occurred
    if callback.condition === nothing
      new_t = zero(typeof(integrator.t))
      min_event_idx = event_idx[1]
    else
      if callback.interp_points!=0
        top_t = ts[interp_index] # Top at the smallest
        bottom_t = ts[prev_sign_index]
      else
        top_t = integrator.t
        bottom_t = integrator.tprev
      end
      if callback.rootfind != NoRootFind && !isdiscrete(integrator.alg)
        min_t = nextfloat(top_t)
        min_event_idx = -1
        for idx in event_idx
          zero_func(abst, p=nothing) = ArrayInterface.allowed_getindex(get_condition(integrator, callback, abst),idx)
          if zero_func(top_t) == 0
            Θ = top_t
          else
            if integrator.event_last_time == counter &&
              integrator.vector_event_last_time == idx &&
              abs(zero_func(bottom_t)) <= 100abs(integrator.last_event_error) &&
              prev_sign_index == 1

              # Determined that there is an event by derivative
              # But floating point error may make the end point negative

              sign_top = sign(zero_func(top_t))
              diff_t = integrator.tdir * 2eps(bottom_t)
              bottom_t += diff_t
              iter = 1
              # This check should match the same check in bisection
              while sign(zero_func(bottom_t)) * sign_top >= zero(sign_top) && iter < 12
                diff_t *= 5
                bottom_t = integrator.tprev + diff_t
                iter += 1
              end
              iter == 12 && error("Double callback crossing floating pointer reducer errored. Report this issue.")
            end
            Θ = bisection(zero_func, (bottom_t, top_t), isone(integrator.tdir), callback.rootfind, callback.abstol, callback.reltol)
            if integrator.tdir * Θ < integrator.tdir * min_t
              integrator.last_event_error = ODE_DEFAULT_NORM(zero_func(Θ), Θ)
            end
          end
          if integrator.tdir * Θ < integrator.tdir * min_t
            min_event_idx = idx
            min_t = Θ
          end
        end
        #Θ = prevfloat(...)
        # prevfloat guerentees that the new time is either 1 floating point
        # numbers just before the event or directly at zero, but not after.
        # If there's a barrier which is never supposed to be crossed,
        # then this will ensure that
        # The item never leaves the domain. Otherwise Roots.jl can return
        # a float which is slightly after, making it out of the domain, causing
        # havoc.
        new_t = min_t -integrator.tprev
      elseif interp_index != callback.interp_points && !isdiscrete(integrator.alg)
        new_t = ts[interp_index] - integrator.tprev
        min_event_idx = event_idx[1]
      else
        # If no solve and no interpolants, just use endpoint
        new_t = integrator.dt
        min_event_idx = event_idx[1]
      end
    end
  else
    new_t = zero(typeof(integrator.t))
    min_event_idx = 1
  end

  new_t,ArrayInterface.allowed_getindex(prev_sign,min_event_idx),event_occurred,min_event_idx
end

function apply_callback!(integrator,callback::Union{ContinuousCallback,VectorContinuousCallback},cb_time,prev_sign,event_idx)

  if isadaptive(integrator)
    set_proposed_dt!(integrator, integrator.tdir * max(nextfloat(integrator.opts.dtmin), integrator.tdir * callback.dtrelax * integrator.dt))
  end

  change_t_via_interpolation!(integrator,integrator.tprev+cb_time)

  # handle saveat
  _, savedexactly = savevalues!(integrator)
  saved_in_cb = true

  @inbounds if callback.save_positions[1]
    # if already saved then skip saving
    savedexactly || savevalues!(integrator,true)
  end

  integrator.u_modified = true

  if prev_sign < 0
    if callback.affect! === nothing
      integrator.u_modified = false
    else
      callback isa VectorContinuousCallback ? callback.affect!(integrator,event_idx) : callback.affect!(integrator)
    end
  elseif prev_sign > 0
    if callback.affect_neg! === nothing
      integrator.u_modified = false
    else
      callback isa VectorContinuousCallback ? callback.affect_neg!(integrator,event_idx) : callback.affect_neg!(integrator)
    end
  end

  if integrator.u_modified
    reeval_internals_due_to_modification!(integrator)
    @inbounds if callback.save_positions[2]
      savevalues!(integrator,true)
      saved_in_cb = true
    end
    return true,saved_in_cb
  end
  false,saved_in_cb
end

#Base Case: Just one
@inline function apply_discrete_callback!(integrator,callback::DiscreteCallback)
  saved_in_cb = false
  if callback.condition(integrator.u,integrator.t,integrator)
    # handle saveat
    _, savedexactly = savevalues!(integrator)
    saved_in_cb = true
    @inbounds if callback.save_positions[1]
      # if already saved then skip saving
      savedexactly || savevalues!(integrator,true)
    end
    integrator.u_modified = true
    callback.affect!(integrator)
    @inbounds if callback.save_positions[2]
      savevalues!(integrator,true)
      saved_in_cb = true
    end
  end
  integrator.sol.destats.ncondition += 1
  integrator.u_modified,saved_in_cb
end

#Starting: Get bool from first and do next
@inline function apply_discrete_callback!(integrator,callback::DiscreteCallback,args...)
  apply_discrete_callback!(integrator,apply_discrete_callback!(integrator,callback)...,args...)
end

@inline function apply_discrete_callback!(integrator,discrete_modified::Bool,saved_in_cb::Bool,callback::DiscreteCallback,args...)
  bool,saved_in_cb2 = apply_discrete_callback!(integrator,apply_discrete_callback!(integrator,callback)...,args...)
  discrete_modified || bool, saved_in_cb || saved_in_cb2
end

@inline function apply_discrete_callback!(integrator,discrete_modified::Bool,saved_in_cb::Bool,callback::DiscreteCallback)
  bool,saved_in_cb2 = apply_discrete_callback!(integrator,callback)
  discrete_modified || bool, saved_in_cb || saved_in_cb2
end

function max_vector_callback_length(cs::CallbackSet)
  continuous_callbacks = cs.continuous_callbacks
  maxlen_cb = nothing
  maxlen = -1
  for cb in continuous_callbacks
    if cb isa VectorContinuousCallback && cb.len > maxlen
      maxlen = cb.len
      maxlen_cb = cb
    end
  end
  maxlen_cb
end

"""
$(TYPEDEF)
"""
mutable struct CallbackCache{conditionType,signType}
  tmp_condition::conditionType
  previous_condition::conditionType
  next_sign::signType
  prev_sign::signType
end

function CallbackCache(u,max_len,::Type{conditionType},::Type{signType}) where {conditionType,signType}
    tmp_condition = similar(u, conditionType, max_len)
    previous_condition = similar(u, conditionType, max_len)
    next_sign = similar(u, signType, max_len)
    prev_sign = similar(u, signType, max_len)
    CallbackCache(tmp_condition,previous_condition,next_sign,prev_sign)
end

function CallbackCache(max_len,::Type{conditionType},::Type{signType}) where {conditionType,signType}
    tmp_condition = zeros(conditionType, max_len)
    previous_condition = zeros(conditionType, max_len)
    next_sign = zeros(signType, max_len)
    prev_sign = zeros(signType, max_len)
    CallbackCache(tmp_condition,previous_condition,next_sign,prev_sign)
end
