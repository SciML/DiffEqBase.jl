# f(t,u,du,res) = 0
"""
$(TYPEDEF)
"""
struct DAEProblem{uType,duType,tType,isinplace,P,F,K,D} <: AbstractDAEProblem{uType,duType,tType,isinplace}
  f::F
  du0::duType
  u0::uType
  tspan::tType
  p::P
  kwargs::K
  differential_vars::D
  @add_kwonly function DAEProblem{iip}(f::AbstractDAEFunction{iip},
                      du0,u0,tspan,p=NullParameters();
                      differential_vars = nothing,
                      kwargs...) where {iip}
    # Defend against external solvers like Sundials breaking on non-uniform input dimensions.
    size(du0) == size(u0) || throw(ArgumentError("Sizes of u0 and du0 must be the same."))
    if !isnothing(differential_vars)
        size(u0) == size(differential_vars) || throw(ArgumentError("Sizes of u0 and differential_vars must be the same."))
    end
    _tspan = promote_tspan(tspan)
    new{typeof(u0),typeof(du0),typeof(_tspan),
               isinplace(f),typeof(p),
               typeof(f),typeof(kwargs),
               typeof(differential_vars)}(
               f,du0,u0,_tspan,p,
               kwargs,differential_vars)
  end

  function DAEProblem{iip}(f,du0,u0,tspan,p=NullParameters();kwargs...) where {iip}
    DAEProblem(convert(DAEFunction{iip},f),du0,u0,tspan,p;kwargs...)
  end
end

function DAEProblem(f::AbstractDAEFunction,du0,u0,tspan,p=NullParameters();kwargs...)
  DAEProblem{isinplace(f)}(f,du0,u0,tspan,p;kwargs...)
end

function DAEProblem(f,du0,u0,tspan,p=NullParameters();kwargs...)
  DAEProblem(convert(DAEFunction,f),du0,u0,tspan,p;kwargs...)
end
