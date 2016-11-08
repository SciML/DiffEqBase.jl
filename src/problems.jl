"""
`ODEProblem`

Wraps the data which defines an ODE problem

```math
\\frac{du}{dt} = f(t,u)
```

with initial condition ``u0``.

### Constructors

`ODEProblem(f,u0,tspan)` : Defines the ODE with the specified functions and
defines the solution if analytic is given.

### Fields

* `f`: The drift function in the ODE.
* `u0`: The initial condition.
* `isinplace`: Determines whether the function `f` uses the in-place syntax `f(t,u,du)`
  or not, `f(t,u)`
* `tspan`: The timespan for the problem.

"""
type ODEProblem{uType,tType,isinplace} <: AbstractODEProblem{uType,tType,isinplace}
  f::Base.Callable
  u0::uType
  tspan::Vector{tType}
end

"""
`ODETestProblem`

Wraps the data which defines an ODE problem

```math
\\frac{du}{dt} = f(t,u)
```

with initial condition ``u0``.

### Constructors

`ODETestProblem(f,u0,analytic,tspan=[0;1])` : Defines the ODE with the specified functions and
defines the solution if analytic is given.

### Fields

* `f`: The drift function in the ODE.
* `u0`: The initial condition.
* `analytic`: A function which describes the solution.
* `numvars`: The number of variables in the system

"""
type ODETestProblem{uType,tType,isinplace} <: AbstractODEProblem{uType,tType,isinplace}
  f::Base.Callable
  u0::uType
  analytic::Base.Callable
  tspan::Vector{tType}
end

function ODEProblem(f::Base.Callable,u0,tspan)
  isinplace = numparameters(f)>=3
  ODEProblem{typeof(u0),eltype(tspan),Val{isinplace}}(f,u0,tspan)
end

function ODETestProblem(f::Base.Callable,u0,analytic,tspan=[0,1.0])
  isinplace = numparameters(f)>=3
  ODETestProblem{typeof(u0),eltype(tspan),Val{isinplace}}(f,u0,analytic,tspan)
end

function print{uType,tType,isinplace}(io::IO, prob::AbstractODEProblem{uType,tType,Val{isinplace}})
  println(io,"AbstractODEProblem")
  println(io,"Independent Variable Type: $uType")
  println(io,"Depdendent Variable Type: $tType")
  println(io,"Function is in-place? $isinplace")
  nothing
end

function show{uType,tType,isinplace}(io::IO,prob::AbstractODEProblem{uType,tType,Val{isinplace}})
  println(io,"AbstractODEProblem{$uType,$tType,$isinplace}")
  nothing
end
