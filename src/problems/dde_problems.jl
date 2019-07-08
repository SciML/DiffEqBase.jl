struct DDEProblem{uType,tType,lType,lType2,isinplace,P,F,H,C} <:
                          AbstractDDEProblem{uType,tType,lType,isinplace}
  f::F
  u0::uType
  h::H
  tspan::tType
  p::P
  constant_lags::lType
  dependent_lags::lType2
  callback::C
  neutral::Bool
  order_discontinuity_t0::Int

  @add_kwonly function DDEProblem{iip}(f::AbstractDDEFunction{iip}, u0, h, tspan, p = nothing;
                                       constant_lags = (),
                                       dependent_lags = (),
                                       neutral = f.mass_matrix !== I && det(f.mass_matrix) != 1,
                                       order_discontinuity_t0 = 0,
                                       callback = nothing) where {iip}
    _tspan = promote_tspan(tspan)
    new{typeof(u0),typeof(_tspan),typeof(constant_lags),typeof(dependent_lags),isinplace(f),
        typeof(p),typeof(f),typeof(h),typeof(callback)}(
          f, u0, h, _tspan, p, constant_lags, dependent_lags, callback, neutral,
          order_discontinuity_t0)
  end

  function DDEProblem{iip}(f::AbstractDDEFunction{iip}, h, tspan::Tuple, p = nothing;
                           order_discontinuity_t0 = 1, kwargs...) where iip
    DDEProblem{iip}(f, h(p, first(tspan)), h, tspan, p;
                    order_discontinuity_t0 = max(1, order_discontinuity_t0), kwargs...)
  end

  function DDEProblem{iip}(f, args...; kwargs...) where iip
    DDEProblem{iip}(convert(DDEFunction{iip}, f), args...; kwargs...)
  end
end

DDEProblem(f, args...; kwargs...) =
  DDEProblem(convert(DDEFunction, f), args...; kwargs...)

DDEProblem(f::AbstractDDEFunction, args...; kwargs...) =
  DDEProblem{isinplace(f)}(f, args...; kwargs...)
