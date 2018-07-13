# Necessary to have initialize set u_modified to false if all don't do anything
# otherwise unnecessary save
INITIALIZE_DEFAULT(cb,u,t,integrator) = u_modified!(integrator, false)

struct ContinuousCallback{F1,F2,F3,F4,T,T2,I} <: AbstractContinuousCallback
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
                   abstol=10eps(),reltol=0) = ContinuousCallback(
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
                   abstol=10eps(),reltol=0)

 ContinuousCallback(
            condition,affect!,affect_neg!,initialize,idxs,
            rootfind,interp_points,
            save_positions,abstol,reltol)

end

struct DiscreteCallback{F1,F2,F3} <: AbstractDiscreteCallback
  condition::F1
  affect!::F2
  initialize::F3
  save_positions::Tuple{Bool,Bool}
end
DiscreteCallback(condition,affect!;
        initialize = INITIALIZE_DEFAULT,save_positions=(true,true)) = DiscreteCallback(condition,affect!,initialize,save_positions)

# DiscreteCallback(condition,affect!,save_positions) = DiscreteCallback(condition,affect!,save_positions)

struct CallbackSet{T1<:Tuple,T2<:Tuple} <: DECallback
  continuous_callbacks::T1
  discrete_callbacks::T2
end

CallbackSet(callback::AbstractDiscreteCallback) = CallbackSet((),(callback,))
CallbackSet(callback::AbstractContinuousCallback) = CallbackSet((callback,),())
CallbackSet() = CallbackSet((),())
CallbackSet(cb::Void) = CallbackSet()

# For Varargs, use recursion to make it type-stable
CallbackSet(callbacks::Union{DECallback,Void}...) = CallbackSet(split_callbacks((), (), callbacks...)...)

@inline split_callbacks(cs, ds) = cs, ds
@inline split_callbacks(cs, ds, c::Void, args...) = split_callbacks(cs, ds, args...)
@inline split_callbacks(cs, ds, c::AbstractContinuousCallback, args...) = split_callbacks((cs..., c), ds, args...)
@inline split_callbacks(cs, ds, d::AbstractDiscreteCallback, args...) = split_callbacks(cs, (ds..., d), args...)
@inline function split_callbacks(cs, ds, d::CallbackSet, args...)
  split_callbacks((cs...,d.continuous_callbacks...), (ds..., d.discrete_callbacks...), args...)
end

# Recursively apply initialize! and return whether any modified u
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

# Helpers
Base.isempty(cb::CallbackSet) = isempty(cb.continuous_callbacks) && isempty(cb.discrete_callbacks)
Base.isempty(cb::AbstractContinuousCallback) = false
Base.isempty(cb::AbstractDiscreteCallback) = false
