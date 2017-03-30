INITIALIZE_DEFAULT(cb,t,u,integrator) = nothing

immutable ContinuousCallback{F1,F2,F3,F4,T,T2,I} <: DECallback
  condition::F1
  affect!::F2
  affect_neg!::F3
  initialize::F4
  idxs::I
  rootfind::Bool
  interp_points::Int
  save_positions::Tuple{Bool,Bool}
  abstol::T
  reltol::T2
end

ContinuousCallback(condition,affect!,affect_neg!;
                   initialize = INITIALIZE_DEFAULT,
                   idxs = nothing,
                   rootfind=true,
                   save_positions=(true,true),
                   interp_points=10,
                   abstol=1e-12,reltol=0) = ContinuousCallback(
                              condition,affect!,affect_neg!,initialize,
                              idxs,
                              rootfind,interp_points,
                              save_positions,abstol,reltol)

function ContinuousCallback(condition,affect!;
                   initialize = INITIALIZE_DEFAULT,
                   idxs = nothing,
                   rootfind=true,
                   save_positions=(true,true),
                   affect_neg! = affect!,
                   interp_points=10,
                   abstol=1e-12,reltol=0)

 ContinuousCallback(
            condition,affect!,affect_neg!,initialize,idxs,
            rootfind,interp_points,
            save_positions,abstol,reltol)

end

immutable DiscreteCallback{F1,F2,F3} <: DECallback
  condition::F1
  affect!::F2
  initialize::F3
  save_positions::Tuple{Bool,Bool}
end
DiscreteCallback(condition,affect!;
        initialize = INITIALIZE_DEFAULT,save_positions=(true,true)) = DiscreteCallback(condition,affect!,initialize,save_positions)

# DiscreteCallback(condition,affect!,save_positions) = DiscreteCallback(condition,affect!,save_positions)

immutable CallbackSet{T1,T2} <: DECallback
  continuous_callbacks::T1
  discrete_callbacks::T2
end

CallbackSet(callback::DiscreteCallback) = CallbackSet((),(callback,))
CallbackSet(callback::ContinuousCallback) = CallbackSet((callback,),())
CallbackSet() = CallbackSet((),())
CallbackSet(cb::Void) = CallbackSet()

# For Varargs, use recursion to make it type-stable

CallbackSet(callbacks::DECallback...) = CallbackSet(split_callbacks((), (), callbacks...)...)

@inline split_callbacks(cs, ds) = cs, ds
@inline split_callbacks(cs, ds, c::ContinuousCallback, args...) = split_callbacks((cs..., c), ds, args...)
@inline split_callbacks(cs, ds, d::DiscreteCallback, args...) = split_callbacks(cs, (ds..., d), args...)
@inline function split_callbacks(cs, ds, d::CallbackSet, args...)
  split_callbacks((cs...,d.continuous_callbacks...), (ds..., d.discrete_callbacks...), args...)
end

initialize!(cb::CallbackSet,t,u,integrator::DEIntegrator) = initialize!(t,u,integrator,cb.continuous_callbacks...,cb.discrete_callbacks...)
initialize!{T}(cb::CallbackSet{Tuple{},T},t,u,integrator::DEIntegrator) = initialize!(t,u,integrator,cb.discrete_callbacks...)
initialize!{T}(cb::CallbackSet{T,Tuple{}},t,u,integrator::DEIntegrator) = initialize!(t,u,integrator,cb.continuous_callbacks...)
initialize!(cb::CallbackSet{Tuple{},Tuple{}},t,u,integrator::DEIntegrator) = nothing
initialize!(t,u,integrator::DEIntegrator,c::DECallback,cs::DECallback...) = (c.initialize(c,t,u,integrator); initialize!(t,u,integrator,cs...))
initialize!(t,u,integrator::DEIntegrator,c::DECallback) = c.initialize(c,t,u,integrator)
