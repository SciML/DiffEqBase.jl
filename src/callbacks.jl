"""
    initialize!(cb::CallbackSet,u,t,integrator::DEIntegrator)

Recursively apply `initialize!` and return whether any modified u
"""
function initialize!(cb::CallbackSet, u, t, integrator::DEIntegrator)
    initialize!(u, t, integrator, false, cb.continuous_callbacks...,
        cb.discrete_callbacks...)
end
initialize!(cb::CallbackSet{Tuple{}, Tuple{}}, u, t, integrator::DEIntegrator) = false
function initialize!(u, t, integrator::DEIntegrator, any_modified::Bool,
        c::DECallback, cs::DECallback...)
    c.initialize(c, u, t, integrator)
    initialize!(u, t, integrator, any_modified || integrator.u_modified, cs...)
end
function initialize!(u, t, integrator::DEIntegrator, any_modified::Bool,
        c::DECallback)
    c.initialize(c, u, t, integrator)
    any_modified || integrator.u_modified
end

"""
    finalize!(cb::CallbackSet,u,t,integrator::DEIntegrator)

Recursively apply `finalize!` and return whether any modified u
"""
function finalize!(cb::CallbackSet, u, t, integrator::DEIntegrator)
    finalize!(u, t, integrator, false, cb.continuous_callbacks..., cb.discrete_callbacks...)
end
finalize!(cb::CallbackSet{Tuple{}, Tuple{}}, u, t, integrator::DEIntegrator) = false
function finalize!(u, t, integrator::DEIntegrator, any_modified::Bool,
        c::DECallback, cs::DECallback...)
    c.finalize(c, u, t, integrator)
    finalize!(u, t, integrator, any_modified || integrator.u_modified, cs...)
end
function finalize!(u, t, integrator::DEIntegrator, any_modified::Bool,
        c::DECallback)
    c.finalize(c, u, t, integrator)
    any_modified || integrator.u_modified
end

# Helpers
function Base.isempty(cb::CallbackSet)
    isempty(cb.continuous_callbacks) && isempty(cb.discrete_callbacks)
end
Base.isempty(cb::AbstractContinuousCallback) = false
Base.isempty(cb::AbstractDiscreteCallback) = false

has_continuous_callback(cb::DiscreteCallback) = false
has_continuous_callback(cb::ContinuousCallback) = true
has_continuous_callback(cb::VectorContinuousCallback) = true
has_continuous_callback(cb::CallbackSet) = !isempty(cb.continuous_callbacks)
has_continuous_callback(cb::Nothing) = false

# Callback handling

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
                integrator(tmp, abst, Val{0})
            else
                integrator(tmp, abst, Val{0}, idxs = callback.idxs)
            end
        else
            if callback.idxs === nothing
                tmp = integrator(abst, Val{0})
            else
                tmp = integrator(abst, Val{0}, idxs = callback.idxs)
            end
        end
        # ismutable && !(callback.idxs isa Number) ? integrator(tmp,abst,Val{0},idxs=callback.idxs) :
        #                                                 tmp = integrator(abst,Val{0},idxs=callback.idxs)
    end
    integrator.sol.stats.ncondition += 1
    if callback isa VectorContinuousCallback
        callback.condition(
            @view(integrator.callback_cache.tmp_condition[1:(callback.len)]),
            tmp, abst, integrator)
        return @view(integrator.callback_cache.tmp_condition[1:(callback.len)])
    else
        return callback.condition(tmp, abst, integrator)
    end
end

# Use a generated function for type stability even when many callbacks are given
@inline function find_first_continuous_callback(integrator,
        callbacks::Vararg{
            AbstractContinuousCallback,
            N}) where {N}
    find_first_continuous_callback(integrator, tuple(callbacks...))
end
@generated function find_first_continuous_callback(integrator,
        callbacks::NTuple{N,
            AbstractContinuousCallback
        }) where {N}
    ex = quote
        tmin, upcrossing,
        event_occurred, event_idx = find_callback_time(integrator,
            callbacks[1], 1)
        identified_idx = 1
    end
    for i in 2:N
        ex = quote
            $ex
            tmin2, upcrossing2,
            event_occurred2, event_idx2 = find_callback_time(
                integrator,
                callbacks[$i],
                $i)
            if event_occurred2 && (tmin2 < tmin || !event_occurred)
                tmin = tmin2
                upcrossing = upcrossing2
                event_occurred = true
                event_idx = event_idx2
                identified_idx = $i
            end
        end
    end
    ex = quote
        $ex
        return tmin, upcrossing, event_occurred, event_idx, identified_idx, $N
    end
    ex
end

@inline function determine_event_occurrence(integrator, callback::VectorContinuousCallback,
        counter)
    event_occurred = false
    if callback.interp_points != 0
        addsteps!(integrator)
    end

    ts = range(integrator.tprev, stop = integrator.t, length = callback.interp_points)

    #=
    # Faster but can be inaccurate
    if callback.interp_points > 1
      dt = (integrator.t - integrator.tprev) / (callback.interp_points-1)
    else
      dt = integrator.dt
    end
    ts = integrator.tprev:dt:integrator.t
    =#

    interp_index = 0
    # Check if the event occurred
    previous_condition = @views(integrator.callback_cache.previous_condition[1:(callback.len)])

    if callback.idxs === nothing
        callback.condition(previous_condition, integrator.uprev, integrator.tprev,
            integrator)
    else
        callback.condition(previous_condition, integrator.uprev[callback.idxs],
            integrator.tprev, integrator)
    end
    integrator.sol.stats.ncondition += 1

    ivec = integrator.vector_event_last_time
    prev_sign = @view(integrator.callback_cache.prev_sign[1:(callback.len)])
    next_sign = @view(integrator.callback_cache.next_sign[1:(callback.len)])

    if integrator.event_last_time == counter &&
       minimum(ODE_DEFAULT_NORM(
        ArrayInterface.allowed_getindex(previous_condition,
            ivec), integrator.t)) <=
       100ODE_DEFAULT_NORM(integrator.last_event_error, integrator.t)

        # If there was a previous event, utilize the derivative at the start to
        # chose the previous sign. If the derivative is positive at tprev, then
        # we treat `prev_sign` as negative, and if the derivative is negative then we
        # treat `prev_sign` as positive, regardless of the positivity/negativity
        # of the true value due to it being =0 sans floating point issues.

        # Only due this if the discontinuity did not move it far away from an event
        # Since near even we use direction instead of location to reset

        if callback.interp_points == 0
            addsteps!(integrator)
        end

        # Evaluate condition slightly in future
        abst = integrator.tprev + integrator.dt * callback.repeat_nudge
        tmp_condition = get_condition(integrator, callback, abst)
        @. prev_sign = sign(previous_condition)
        prev_sign[ivec] = sign(tmp_condition[ivec])
    else
        @. prev_sign = sign(previous_condition)
    end

    prev_sign_index = 1
    abst = integrator.t
    next_condition = get_condition(integrator, callback, abst)
    @. next_sign = sign(next_condition)

    event_idx = findall_events!(next_sign, callback.affect!, callback.affect_neg!,
        prev_sign)
    if sum(event_idx) != 0
        event_occurred = true
        interp_index = callback.interp_points
    end

    if callback.interp_points != 0 && !isdiscrete(integrator.alg) &&
       sum(event_idx) != length(event_idx) # Use the interpolants for safety checking
        fallback = true
        for i in 2:length(ts)
            abst = ts[i]
            copyto!(next_sign, get_condition(integrator, callback, abst))
            _event_idx = findall_events!(next_sign, callback.affect!, callback.affect_neg!,
                prev_sign)
            if sum(_event_idx) != 0
                event_occurred = true
                event_idx = _event_idx
                interp_index = i
                fallback = false
                break
            else
                prev_sign_index = i
            end
        end

        if fallback
            # If you get here, then you need to reset the event_idx to the
            # non-interpolated version

            abst = integrator.t
            next_condition = get_condition(integrator, callback, abst)
            @. next_sign = sign(next_condition)
            event_idx = findall_events!(next_sign, callback.affect!, callback.affect_neg!,
                prev_sign)
            interp_index = callback.interp_points
        end
    end

    event_occurred, interp_index, ts, prev_sign, prev_sign_index, event_idx
end

@inline function determine_event_occurrence(integrator, callback::ContinuousCallback,
        counter)
    event_occurred = false
    if callback.interp_points != 0
        addsteps!(integrator)
    end

    ts = range(integrator.tprev, stop = integrator.t, length = callback.interp_points)

    #=
    # Faster but can be inaccurate
    if callback.interp_points > 1
      dt = (integrator.t - integrator.tprev) / (callback.interp_points-1)
    else
      dt = integrator.dt
    end
    ts = integrator.tprev:dt:integrator.t
    =#

    interp_index = 0

    # Check if the event occurred
    if callback.idxs === nothing
        previous_condition = callback.condition(integrator.uprev, integrator.tprev,
            integrator)
    else
        @views previous_condition = callback.condition(integrator.uprev[callback.idxs],
            integrator.tprev, integrator)
    end
    integrator.sol.stats.ncondition += 1

    prev_sign = 0.0
    next_sign = 0.0

    if integrator.event_last_time == counter &&
       minimum(ODE_DEFAULT_NORM(previous_condition, integrator.t)) <=
       100ODE_DEFAULT_NORM(integrator.last_event_error, integrator.t)

        # If there was a previous event, utilize the derivative at the start to
        # chose the previous sign. If the derivative is positive at tprev, then
        # we treat `prev_sign` as negative, and if the derivative is negative then we
        # treat `prev_sign` as positive, regardless of the positivity/negativity
        # of the true value due to it being =0 sans floating point issues.

        # Only due this if the discontinuity did not move it far away from an event
        # Since near even we use direction instead of location to reset

        if callback.interp_points == 0
            addsteps!(integrator)
        end

        # Evaluate condition slightly in future
        abst = integrator.tprev + integrator.dt * callback.repeat_nudge
        tmp_condition = get_condition(integrator, callback, abst)
        prev_sign = sign(tmp_condition)
    else
        prev_sign = sign(previous_condition)
    end

    prev_sign_index = 1
    abst = integrator.t
    next_condition = get_condition(integrator, callback, abst)
    next_sign = sign(next_condition)

    if ((prev_sign < 0 && callback.affect! !== nothing) ||
        (prev_sign > 0 && callback.affect_neg! !== nothing)) && prev_sign * next_sign <= 0
        event_occurred = true
        interp_index = callback.interp_points
    elseif callback.interp_points != 0 && !isdiscrete(integrator.alg) # Use the interpolants for safety checking
        for i in 2:length(ts)
            abst = ts[i]
            new_sign = get_condition(integrator, callback, abst)
            if ((prev_sign < 0 && callback.affect! !== nothing) ||
                (prev_sign > 0 && callback.affect_neg! !== nothing)) &&
               prev_sign * new_sign < 0
                event_occurred = true
                interp_index = i
                break
            else
                prev_sign_index = i
            end
        end
    end
    event_idx = 1

    event_occurred, interp_index, ts, prev_sign, prev_sign_index, event_idx
end

# rough implementation, needs multiple type handling
# always ensures that if r = bisection(f, (x0, x1))
# then either f(nextfloat(r)) == 0 or f(nextfloat(r)) * f(r) < 0
# note: not really using bisection - uses the ITP method
function bisection(
        f, tup, t_forward::Bool, rootfind::SciMLBase.RootfindOpt, abstol, reltol;
        maxiters = 1000)
    if rootfind == SciMLBase.LeftRootFind
        solve(IntervalNonlinearProblem{false}(f, tup),
            InternalITP(), abstol = abstol,
            reltol = reltol).left
    else
        solve(IntervalNonlinearProblem{false}(f, tup),
            InternalITP(), abstol = abstol,
            reltol = reltol).right
    end
end

"""
findall_events!(next_sign,affect!,affect_neg!,prev_sign)

Modifies `next_sign` to be an array of booleans for if there is a sign change
in the interval between prev_sign and next_sign
"""
function findall_events!(next_sign::Union{Array, SubArray}, affect!::F1, affect_neg!::F2,
        prev_sign::Union{Array, SubArray}) where {F1, F2}
    @inbounds for i in 1:length(prev_sign)
        next_sign[i] = ((prev_sign[i] < 0 && affect! !== nothing) ||
                        (prev_sign[i] > 0 && affect_neg! !== nothing)) &&
                       prev_sign[i] * next_sign[i] <= 0
    end
    next_sign
end

function findall_events!(next_sign, affect!::F1, affect_neg!::F2, prev_sign) where {F1, F2}
    hasaffect::Bool = affect! !== nothing
    hasaffectneg::Bool = affect_neg! !== nothing
    f = (n, p) -> ((p < 0 && hasaffect) || (p > 0 && hasaffectneg)) && p * n <= 0
    A = map!(f, next_sign, next_sign, prev_sign)
    next_sign
end

function find_callback_time(integrator, callback::ContinuousCallback, counter)
    event_occurred, interp_index, ts, prev_sign,
    prev_sign_index, event_idx = determine_event_occurrence(
        integrator,
        callback,
        counter)
    if event_occurred
        if callback.condition === nothing
            new_t = zero(typeof(integrator.t))
        else
            if callback.interp_points != 0
                top_t = ts[interp_index] # Top at the smallest
                bottom_t = ts[prev_sign_index]
            else
                top_t = integrator.t
                bottom_t = integrator.tprev
            end
            if callback.rootfind != SciMLBase.NoRootFind && !isdiscrete(integrator.alg)
                zero_func(abst, p = nothing) = get_condition(integrator, callback, abst)
                if zero_func(top_t) == 0
                    Θ = top_t
                else
                    if integrator.event_last_time == counter &&
                       abs(zero_func(bottom_t)) <= 100abs(integrator.last_event_error) &&
                       prev_sign_index == 1

                        # Determined that there is an event by derivative
                        # But floating point error may make the end point negative

                        bottom_t += integrator.dt * callback.repeat_nudge
                        sign_top = sign(zero_func(top_t))
                        sign(zero_func(bottom_t)) * sign_top >= zero(sign_top) &&
                            error("Double callback crossing floating pointer reducer errored. Report this issue.")
                    end
                    Θ = bisection(zero_func, (bottom_t, top_t), isone(integrator.tdir),
                        callback.rootfind, callback.abstol, callback.reltol)
                    integrator.last_event_error = DiffEqBase.value(ODE_DEFAULT_NORM(
                        zero_func(Θ), Θ))
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

    new_t, prev_sign, event_occurred, event_idx
end

function find_callback_time(integrator, callback::VectorContinuousCallback, counter)
    event_occurred, interp_index, ts, prev_sign,
    prev_sign_index, event_idx = determine_event_occurrence(
        integrator,
        callback,
        counter)
    if event_occurred
        if callback.condition === nothing
            new_t = zero(typeof(integrator.t))
            min_event_idx = findfirst(isequal(1), event_idx)
        else
            if callback.interp_points != 0
                top_t = ts[interp_index] # Top at the smallest
                bottom_t = ts[prev_sign_index]
            else
                top_t = integrator.t
                bottom_t = integrator.tprev
            end
            if callback.rootfind != SciMLBase.NoRootFind && !isdiscrete(integrator.alg)
                min_t = isone(integrator.tdir) ? nextfloat(top_t) : prevfloat(top_t)
                min_event_idx = -1
                for idx in 1:length(event_idx)
                    if ArrayInterface.allowed_getindex(event_idx, idx) != 0
                        function zero_func(abst, p = nothing)
                            ArrayInterface.allowed_getindex(
                                get_condition(integrator,
                                    callback,
                                    abst), idx)
                        end
                        if zero_func(top_t) == 0
                            Θ = top_t
                        else
                            if integrator.event_last_time == counter &&
                               integrator.vector_event_last_time == idx &&
                               abs(zero_func(bottom_t)) <=
                               100abs(integrator.last_event_error) &&
                               prev_sign_index == 1

                                # Determined that there is an event by derivative
                                # But floating point error may make the end point negative

                                bottom_t += integrator.dt * callback.repeat_nudge
                                sign_top = sign(zero_func(top_t))
                                sign(zero_func(bottom_t)) * sign_top >= zero(sign_top) &&
                                    error("Double callback crossing floating pointer reducer errored. Report this issue.")
                            end

                            Θ = bisection(zero_func, (bottom_t, top_t),
                                isone(integrator.tdir), callback.rootfind,
                                callback.abstol, callback.reltol)
                            if integrator.tdir * Θ < integrator.tdir * min_t
                                integrator.last_event_error = DiffEqBase.value(ODE_DEFAULT_NORM(
                                    zero_func(Θ), Θ))
                            end
                        end
                        if integrator.tdir * Θ < integrator.tdir * min_t
                            min_event_idx = idx
                            min_t = Θ
                        end
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
                new_t = min_t - integrator.tprev
            elseif interp_index != callback.interp_points && !isdiscrete(integrator.alg)
                new_t = ts[interp_index] - integrator.tprev
                min_event_idx = findfirst(isequal(1), event_idx)
            else
                # If no solve and no interpolants, just use endpoint
                new_t = integrator.dt
                min_event_idx = findfirst(isequal(1), event_idx)
            end
        end
    else
        new_t = zero(typeof(integrator.t))
        min_event_idx = 1
    end

    if event_occurred && min_event_idx < 0
        error("Callback handling failed. Please file an issue with code to reproduce.")
    end

    new_t, ArrayInterface.allowed_getindex(prev_sign, min_event_idx),
    event_occurred::Bool, min_event_idx::Int
end

function apply_callback!(integrator,
        callback::Union{ContinuousCallback, VectorContinuousCallback},
        cb_time, prev_sign, event_idx)
    if isadaptive(integrator)
        set_proposed_dt!(integrator,
            integrator.tdir * max(nextfloat(integrator.opts.dtmin),
                integrator.tdir * callback.dtrelax * integrator.dt))
    end

    change_t_via_interpolation!(
        integrator, integrator.tprev + cb_time, Val{:false}, callback.initializealg)

    # handle saveat
    _, savedexactly = savevalues!(integrator)
    saved_in_cb = true

    @inbounds if callback.save_positions[1]
        # if already saved then skip saving
        savedexactly || savevalues!(integrator, true)
    end

    integrator.u_modified = true

    if prev_sign < 0
        if callback.affect! === nothing
            integrator.u_modified = false
        else
            callback isa VectorContinuousCallback ?
            callback.affect!(integrator, event_idx) : callback.affect!(integrator)
        end
    elseif prev_sign > 0
        if callback.affect_neg! === nothing
            integrator.u_modified = false
        else
            callback isa VectorContinuousCallback ?
            callback.affect_neg!(integrator, event_idx) : callback.affect_neg!(integrator)
        end
    end

    if integrator.u_modified
        reeval_internals_due_to_modification!(
            integrator, callback_initializealg = callback.initializealg)

        @inbounds if callback.save_positions[2]
            savevalues!(integrator, true)
            if !isdefined(integrator.opts, :save_discretes) || integrator.opts.save_discretes
                if callback isa VectorContinuousCallback
                    SciMLBase.save_discretes!(integrator, callback, event_idx)
                else
                    SciMLBase.save_discretes!(integrator, callback)
                end
            end
            saved_in_cb = true
        end
        return true, saved_in_cb
    end
    false, saved_in_cb
end

#Base Case: Just one
@inline function apply_discrete_callback!(integrator, callback::DiscreteCallback)
    saved_in_cb = false
    if callback.condition(integrator.u, integrator.t, integrator)
        # handle saveat
        _, savedexactly = savevalues!(integrator)
        saved_in_cb = true
        @inbounds if callback.save_positions[1]
            # if already saved then skip saving
            savedexactly || savevalues!(integrator, true)
        end
        integrator.u_modified = true
        callback.affect!(integrator)
        if integrator.u_modified
            reeval_internals_due_to_modification!(
                integrator, false, callback_initializealg = callback.initializealg)
        end
        @inbounds if callback.save_positions[2]
            savevalues!(integrator, true)
            if !isdefined(integrator.opts, :save_discretes) || integrator.opts.save_discretes
                SciMLBase.save_discretes!(integrator, callback)
            end
            saved_in_cb = true
        end
    end
    integrator.sol.stats.ncondition += 1
    integrator.u_modified, saved_in_cb
end

#Starting: Get bool from first and do next
@inline function apply_discrete_callback!(integrator, callback::DiscreteCallback, args...)
    apply_discrete_callback!(integrator, apply_discrete_callback!(integrator, callback)...,
        args...)
end

@inline function apply_discrete_callback!(integrator, discrete_modified::Bool,
        saved_in_cb::Bool, callback::DiscreteCallback,
        args...)
    bool,
    saved_in_cb2 = apply_discrete_callback!(integrator,
        apply_discrete_callback!(integrator,
            callback)...,
        args...)
    discrete_modified || bool, saved_in_cb || saved_in_cb2
end

@inline function apply_discrete_callback!(integrator, discrete_modified::Bool,
        saved_in_cb::Bool, callback::DiscreteCallback)
    bool, saved_in_cb2 = apply_discrete_callback!(integrator, callback)
    discrete_modified || bool, saved_in_cb || saved_in_cb2
end

function max_vector_callback_length_int(cs::CallbackSet)
    max_vector_callback_length_int(cs.continuous_callbacks...)
end
max_vector_callback_length_int() = nothing
function max_vector_callback_length_int(continuous_callbacks...)
    all(cb -> cb isa ContinuousCallback, continuous_callbacks) && return nothing
    maxlen = -1
    for cb in continuous_callbacks
        if cb isa VectorContinuousCallback && cb.len > maxlen
            maxlen = cb.len
        end
    end
    maxlen
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
mutable struct CallbackCache{conditionType, signType}
    tmp_condition::conditionType
    previous_condition::conditionType
    next_sign::signType
    prev_sign::signType
end

function CallbackCache(u, max_len, ::Type{conditionType},
        ::Type{signType}) where {conditionType, signType}
    tmp_condition = similar(u, conditionType, max_len)
    previous_condition = similar(u, conditionType, max_len)
    next_sign = similar(u, signType, max_len)
    prev_sign = similar(u, signType, max_len)
    CallbackCache(tmp_condition, previous_condition, next_sign, prev_sign)
end

function CallbackCache(max_len, ::Type{conditionType},
        ::Type{signType}) where {conditionType, signType}
    tmp_condition = zeros(conditionType, max_len)
    previous_condition = zeros(conditionType, max_len)
    next_sign = zeros(signType, max_len)
    prev_sign = zeros(signType, max_len)
    CallbackCache(tmp_condition, previous_condition, next_sign, prev_sign)
end
