"""
Defines a ordinary differential equation (ODE) in explicit form
du/dt = f(t,u).

Input:
- `f(t,u,out)` or `f(t,u)`
- `u0` initial conditions (IC) on `u`
- `tspan` tuple of start and end-time `(t0,tend)`

Extra function, such as the Jacobian can be provided by function overloading:
- `f(::Val{:jac}, t, u, out)` where `jac = df/du + alpha* df/d(du)`

Notes:
- the types of inputs determines the types used in the numerical
  solver.  For instance, using `Int`s in `tspan` will not work.
"""
type ODEProblem{uType,tType,isinplace,F} <: AbstractODEProblem{uType,tType,isinplace,F}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
end

type ODETestProblem{uType,tType,isinplace,F} <: AbstractODETestProblem{uType,tType,isinplace,F}
  f::F
  u0::uType
  analytic::Base.Callable
  tspan::Tuple{tType,tType}
end

function ODEProblem(f,u0,tspan)
  isinplace = numparameters(f)>=3
  ODEProblem{typeof(u0),eltype(tspan),isinplace,typeof(f)}(f,u0,tspan)
end

function ODETestProblem(f,u0,analytic,tspan=(0.0,1.0))
  isinplace = numparameters(f)>=3
  ODETestProblem{typeof(u0),eltype(tspan),isinplace,typeof(f)}(f,u0,analytic,tspan)
end

function print{uType,tType,isinplace,F}(io::IO, prob::AbstractODEProblem{uType,tType,isinplace,F})
  println(io,"AbstractODEProblem")
  println(io,"Independent Variable Type: $uType")
  println(io,"Depdendent Variable Type: $tType")
  println(io,"Function is in-place? $isinplace")
  nothing
end

#=
function show{uType,tType,isinplace,F}(io::IO,prob::AbstractODEProblem{uType,tType,isinplace,F})
  println(io,"AbstractODEProblem{$uType,$tType,$isinplace}")
  nothing
end
=#
