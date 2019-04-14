"""
$(TYPEDEF)

TODO
"""
struct StandardODEProblem end

# Mu' = f
"""
$(TYPEDEF)

Defines an ODE problem.

# Fields

$(FIELDS)
"""
struct ODEProblem{uType,tType,isinplace,P,F,C,PT} <:
               AbstractODEProblem{uType,tType,isinplace}
  """The function in the ODE."""
  f::F
  """The initial condition."""
  u0::uType
  """The timespan for the problem."""
  tspan::tType
  """The parameter values of the ODE function."""
  p::P
  """A callback to be applied to every solver which uses the problem."""
  callback::C
  """TODO"""
  problem_type::PT
  @add_kwonly function ODEProblem{iip}(f::AbstractODEFunction{iip},
                                       u0,tspan,p=nothing,
                                       problem_type=StandardODEProblem();
                                       callback=nothing) where {iip}
    _tspan = promote_tspan(tspan)
    new{typeof(u0),typeof(_tspan),
       isinplace(f),typeof(p),typeof(f),
       typeof(callback),
       typeof(problem_type)}(
       f,u0,_tspan,p,callback,problem_type)
  end

  """
      ODEProblem{isinplace}(f,u0,tspan,p=nothing,callback=CallbackSet())

  Define an ODE problem with the specified function.
  `isinplace` optionally sets whether the function is inplace or not.
  This is determined automatically, but not inferred.
  """
  function ODEProblem{iip}(f,u0,tspan,p=nothing;kwargs...) where {iip}
    ODEProblem(convert(ODEFunction{iip},f),u0,tspan,p;kwargs...)
  end

  function ODEProblem{iip,recompile}(f,u0,tspan,p=nothing;kwargs...) where {iip,recompile}
    if !recompile
      if iip
        ODEProblem{iip}(wrapfun_iip(f,(u0,u0,p,tspan[1])),u0,tspan,p;kwargs...)
      else
        ODEProblem{iip}(wrapfun_oop(f,(u0,p,tspan[1])),u0,tspan,p;kwargs...)
      end
    else
      ODEProblem{iip}(f,u0,tspan,p;kwargs...)
    end
  end
end

"""
    ODEProblem(f::ODEFunction,u0,tspan,p=nothing,callback=CallbackSet())

Define an ODE problem from a [`ODEFunction`](@ref).
"""
function ODEProblem(f::AbstractODEFunction,u0,tspan,args...;kwargs...)
  ODEProblem{isinplace(f)}(f,u0,tspan,args...;kwargs...)
end

function ODEProblem(f,u0,tspan,p=nothing;kwargs...)
  ODEProblem(convert(ODEFunction,f),u0,tspan,p;kwargs...)
end

"""
$(TYPEDEF)

TODO
"""
abstract type AbstractDynamicalODEProblem end

"""
$(TYPEDEF)

TODO
"""
struct DynamicalODEProblem{iip} <: AbstractDynamicalODEProblem end
# u' = f1(v)
# v' = f2(t,u)
"""
    DynamicalODEProblem(f::DynamicalODEFunction,v0,u0,tspan,p=nothing,callback=CallbackSet())

Define a dynamical ODE function from a [`DynamicalODEFunction`](@ref).
"""
function DynamicalODEProblem(f::DynamicalODEFunction,du0,u0,tspan,p=nothing;kwargs...)
  ODEProblem(f,(du0,u0),tspan,p;kwargs...)
end
function DynamicalODEProblem(f1,f2,du0,u0,tspan,p=nothing;kwargs...)
  ODEProblem(DynamicalODEFunction(f1,f2),(du0,u0),tspan,p;kwargs...)
end

"""
    DynamicalODEProblem{isinplace}(f1,f2,v0,u0,tspan,p=nothing,callback=CallbackSet())

Define a dynamical ODE problem from the two functions `f1` and `f2`.

# Arguments
* `f1` and `f2`: The functions in the ODE.
* `v0` and `u0`: The initial conditions.
* `tspan`: The timespan for the problem.
* `p`: Parameter values for `f1` and `f2`.
* `callback`: A callback to be applied to every solver which uses the problem. Defaults to nothing.

`isinplace` optionally sets whether the function is inplace or not.
This is determined automatically, but not inferred.
"""
function DynamicalODEProblem{iip}(f1,f2,du0,u0,tspan,p=nothing;kwargs...) where iip
  ODEProblem(DynamicalODEFunction{iip}(f1,f2),(du0,u0),tspan,p;kwargs...)
end

# u'' = f(t,u,du,ddu)
"""
$(TYPEDEF)

TODO
"""
struct SecondOrderODEProblem{iip} <: AbstractDynamicalODEProblem end
function SecondOrderODEProblem(f,du0,u0,tspan,p=nothing;kwargs...)
  iip = isinplace(f,5)
  SecondOrderODEProblem{iip}(f,du0,u0,tspan,p;kwargs...)
end

"""
    SecondOrderODEProblem{isinplace}(f,du0,u0,tspan,p=nothing,callback=CallbackSet())

Define a second order ODE problem with the specified function.

# Arguments
* `f`: The function for the second derivative.
* `du0`: The initial derivative.
* `u0`: The initial condition.
* `tspan`: The timespan for the problem.
* `callback`: A callback to be applied to every solver which uses the problem. Defaults to nothing.
"""
function SecondOrderODEProblem{iip}(f,du0,u0,tspan,p=nothing;kwargs...) where iip
  if iip
    f2 = function (du,v,u,p,t)
      du .= v
    end
  else
    f2 = function (v,u,p,t)
      v
    end
  end
  _u0 = (du0,u0)
  ODEProblem(DynamicalODEFunction{iip}(f,f2),_u0,tspan,p,
                  SecondOrderODEProblem{iip}();kwargs...)
end
function SecondOrderODEProblem(f::DynamicalODEFunction,du0,u0,tspan,p=nothing;kwargs...)
  iip = isinplace(f.f1, 5)
  _u0 = (du0,u0)
  if f.f2.f == nothing
    if iip
      f2 = function (du,v,u,p,t)
        du .= v
      end
    else
      f2 = function (v,u,p,t)
        v
      end
    end
    return ODEProblem(DynamicalODEFunction{iip}(f.f1,f2;mass_matrix=f.mass_matrix,analytic=f.analytic),_u0,tspan,p,
                  SecondOrderODEProblem{iip}();kwargs...)
  else
    return ODEProblem(DynamicalODEFunction{iip}(f.f1,f.f2;mass_matrix=f.mass_matrix,analytic=f.analytic),_u0,tspan,p,
                  SecondOrderODEProblem{iip}();kwargs...)
  end
end

"""
$(TYPEDEF)

TODO
"""
abstract type AbstractSplitODEProblem end

"""
$(TYPEDEF)

TODO
"""
struct SplitODEProblem{iip} <: AbstractSplitODEProblem end
# u' = Au + f
function SplitODEProblem(f1,f2,u0,tspan,p=nothing;kwargs...)
  f = SplitFunction(f1,f2)
  SplitODEProblem(f,u0,tspan,p;kwargs...)
end
"""
    SplitODEProblem{isinplace}(f1,f2,u0,tspan,p=nothing;kwargs...)

Define a split ODE problem from separate functions `f1` and `f2`.

# Arguments
* `f1`, `f2`: The functions in the ODE.
* `u0`: The initial condition.
* `tspan`: The timespan for the problem.
* `p`: The parameters for the problem.
* `callback`: A callback to be applied to every solver which uses the problem. Defaults to nothing.

The `isinplace parameter` can be omitted and will be determined using the
signature of `f2`. Note that both `f1` and `f2` should support the in-place style if
`isinplace` is true or they should both support the out-of-place style if
`isinplace` is false. You cannot mix up the two styles.
"""
function SplitODEProblem{iip}(f1,f2,u0,tspan,p=nothing;kwargs...) where iip
  f = SplitFunction{iip}(f1,f2)
  SplitODEProblem(f,u0,tspan,p;kwargs...)
end

"""
$(SIGNATURES)

Define a split ODE problem from a [`SplitFunction`](@ref).
"""
SplitODEProblem(f::SplitFunction,u0,tspan,p=nothing;kwargs...) =
  SplitODEProblem{isinplace(f)}(f,u0,tspan,p;kwargs...)
function SplitODEProblem{iip}(f::SplitFunction,u0,tspan,p=nothing;kwargs...) where iip
  if f.cache == nothing && iip
    cache = similar(u0)
    f = SplitFunction{iip}(f.f1, f.f2; mass_matrix=f.mass_matrix,
                     _func_cache=cache, analytic=f.analytic)
  end
  ODEProblem(f,u0,tspan,p,SplitODEProblem{iip}();kwargs...)
end
