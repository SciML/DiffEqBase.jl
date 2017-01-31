immutable ContinuousCallback{F1,F2,F3,T,T2} <: DECallback
  condition::F1
  affect!::F2
  affect_neg!::F3
  rootfind::Bool
  interp_points::Int
  save_positions::Tuple{Bool,Bool}
  abstol::T
  reltol::T2
end

ContinuousCallback(condition,affect!,affect_neg!;
                   rootfind=true,
                   save_positions=(true,true),
                   interp_points=10,
                   abstol=1e-12,reltol=0) = ContinuousCallback(
                              condition,affect!,affect_neg!,
                              rootfind,interp_points,
                              save_positions,abstol,reltol)

function ContinuousCallback(condition,affect!;
                   rootfind=true,
                   save_positions=(true,true),
                   affect_neg! = affect!,
                   interp_points=10,
                   abstol=1e-12,reltol=0)

  if affect_neg! != affect!
    Base.depwarn("the `affect_neg!` keyword has been deprecated for a different constructor. Please consult the documentation.",:ContinuousCallback)
  end

 ContinuousCallback(
            condition,affect!,affect_neg!,
            rootfind,interp_points,
            save_positions,abstol,reltol)

end

immutable DiscreteCallback{F1,F2} <: DECallback
  condition::F1
  affect!::F2
  save_positions::Tuple{Bool,Bool}
end
DiscreteCallback(condition,affect!;save_positions=(true,true)) = DiscreteCallback(condition,affect!,save_positions)

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
@inline split_callbacks(cs, ds, d::CallbackSet, args...) = split_callbacks((cs...,d.continuous_callbacks...), (ds..., d.discrete_callbacks...), args...)
