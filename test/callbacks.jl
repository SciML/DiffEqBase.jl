using DiffEqBase, Test

condition = function (u, t, integrator) # Event when event_f(u,t,k) == 0
    t - 2.95
end

affect! = function (integrator)
    integrator.u = integrator.u + 2
end

rootfind = true
save_positions = (true, true)
callback = ContinuousCallback(condition, affect!; save_positions = save_positions)

cbs = CallbackSet(nothing)
@test typeof(cbs.discrete_callbacks) <: Tuple
@test typeof(cbs.continuous_callbacks) <: Tuple
cbs = CallbackSet(callback, nothing)
@test typeof(cbs.discrete_callbacks) <: Tuple
@test typeof(cbs.continuous_callbacks) <: Tuple
cbs = CallbackSet(callback, CallbackSet())
@test typeof(cbs.discrete_callbacks) <: Tuple
@test typeof(cbs.continuous_callbacks) <: Tuple

condition = function (integrator)
    true
end
affect! = function (integrator) end
save_positions = (true, false)
saving_callback = DiscreteCallback(condition, affect!; save_positions = save_positions)

cbs1 = CallbackSet(callback, saving_callback)

@test length(cbs1.discrete_callbacks) == 1
@test length(cbs1.continuous_callbacks) == 1

cbs2 = CallbackSet(callback)
@test length(cbs2.continuous_callbacks) == 1
@test length(cbs2.discrete_callbacks) == 0

cbs3 = CallbackSet(saving_callback)
@test length(cbs3.discrete_callbacks) == 1
@test length(cbs3.continuous_callbacks) == 0

cbs4 = CallbackSet()
@test length(cbs4.discrete_callbacks) == 0
@test length(cbs4.continuous_callbacks) == 0

cbs5 = CallbackSet(cbs1, cbs2)

@test length(cbs5.discrete_callbacks) == 1
@test length(cbs5.continuous_callbacks) == 2
