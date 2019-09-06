struct SDDEProblem{uType,tType,lType,lType2,isinplace,P,NP,F,G,H,C,ND} <:
                          AbstractSDDEProblem{uType,tType,lType,isinplace,ND}
  f::F
  g::G
  u0::uType
  h::H
  tspan::tType
  p::P
  noise::NP
  constant_lags::lType
  dependent_lags::lType2
  callback::C
  noise_rate_prototype::ND
  seed::UInt64
  neutral::Bool
  order_discontinuity_t0::Int
  # TODO @add_kwonly 
  function SDDEProblem{iip}(f::AbstractSDDEFunction{iip}, g, u0, h, tspan, p = nothing;
                                       noise_rate_prototype = nothing, noise= nothing, seed = UInt64(0),
                                       constant_lags = (), dependent_lags = (),
                                       neutral = f.mass_matrix !== I && det(f.mass_matrix) != 1,
                                       order_discontinuity_t0 = 0,
                                       callback = nothing) where {iip}
    _tspan = promote_tspan(tspan)
    new{typeof(u0),typeof(_tspan),typeof(constant_lags),typeof(dependent_lags),isinplace(f),
        typeof(p),typeof(noise),typeof(f),typeof(g),typeof(h),typeof(callback),typeof(noise_rate_prototype)}(
        f, g, u0, h, _tspan, p, noise, constant_lags, dependent_lags, callback, noise_rate_prototype, seed, neutral, order_discontinuity_t0)
  end

  function SDDEProblem{iip}(f::AbstractSDDEFunction{iip}, g , h, tspan, p = nothing;
                           order_discontinuity_t0 = 1, kwargs...) where iip
    SDDEProblem{iip}(f, g, h(p, first(tspan)), h, tspan, p;
                    order_discontinuity_t0 = max(1, order_discontinuity_t0), kwargs...)
  end

  function SDDEProblem{iip}(f,g, args...; kwargs...) where iip
    SDDEProblem{iip}(convert(SDDEFunction{iip}, f,g),g, args...; kwargs...)
  end
end

SDDEProblem(f, g, args...; kwargs...) =
  SDDEProblem(convert(SDDEFunction, f, g), g, args...; kwargs...)

SDDEProblem(f::AbstractSDDEFunction, args...; kwargs...) =
  SDDEProblem{isinplace(f)}(f, args...; kwargs...)
