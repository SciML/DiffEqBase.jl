# Function Interface

```@meta
CurrentModule = DiffEqBase
```

Functions define the dynamics of differential equation problems. DiffEqBase extends the function types from SciMLBase with additional capabilities.

## Function Types

The following abstract function types are imported from SciMLBase:

- `AbstractODEFunction`: Functions for ODEs
- `AbstractSDEFunction`: Functions for SDEs
- `AbstractRODEFunction`: Functions for RODEs
- `AbstractDDEFunction`: Functions for DDEs
- `AbstractSDDEFunction`: Functions for SDDEs
- `AbstractDAEFunction`: Functions for DAEs
- `AbstractNonlinearFunction`: Functions for nonlinear problems
- `AbstractDiffEqFunction`: Base type for all DE functions

## Function Interface

### In-Place vs Out-of-Place

Functions can be defined in two styles:

**In-place** (mutating):
```julia
function f!(du, u, p, t)
    du[1] = p[1] * u[1]
    du[2] = -p[2] * u[2]
end
```

**Out-of-place** (non-mutating):
```julia
function f(u, p, t)
    [p[1] * u[1], -p[2] * u[2]]
end
```

Use `isinplace(f)` to check the style.

### Creating Function Objects

For advanced features, create function objects:

```julia
# ODE with Jacobian
function f!(du, u, p, t)
    du[1] = p[1] * u[1]
end

function jac!(J, u, p, t)
    J[1,1] = p[1]
end

odefun = ODEFunction(f!; jac=jac!)
```

### Available Properties

Function objects can include:

- `f`: The main dynamics function
- `jac`: Jacobian matrix function
- `tgrad`: Time gradient function
- `paramjac`: Parameter Jacobian function
- `mass_matrix`: Mass matrix for DAEs
- `analytic`: Analytical solution function
- `colorvec`: Sparsity/coloring information

### Checking Properties

- `has_jac(f)`: Check if Jacobian is provided
- `has_tgrad(f)`: Check if time gradient is provided
- `has_paramjac(f)`: Check if parameter Jacobian is provided
- `has_analytic(f)`: Check if analytical solution is provided
- `has_colorvec(f)`: Check if coloring vector is provided

## Specialized Functions

### SDE Functions

For stochastic equations, also define noise:

```julia
function f!(du, u, p, t)
    du[1] = p[1] * u[1]
end

function g!(du, u, p, t)
    du[1] = p[2] * u[1]  # Noise term
end

sdefun = SDEFunction(f!, g!)
```

### DDE Functions

For delay equations, access history:

```julia
function f!(du, u, h, p, t)
    history = h(p, t - p[3])  # Access delayed value
    du[1] = p[1] * (history[1] - u[1])
end
```

### DAE Functions

For differential-algebraic equations:

```julia
function f!(resid, du, u, p, t)
    resid[1] = du[1] - p[1] * u[1]
    resid[2] = u[1] + u[2] - p[2]  # Algebraic constraint
end

daefun = DAEFunction(f!)
```

## Wrappers

DiffEqBase provides wrapper types for function components:

- `TimeGradientWrapper`: Wraps time gradient functions
- `TimeDerivativeWrapper`: Wraps time derivative functions
- `UDerivativeWrapper`: Wraps state derivative functions
- `UJacobianWrapper`: Wraps Jacobian functions
- `ParamJacobianWrapper`: Wraps parameter Jacobian functions
- `JacobianWrapper`: General Jacobian wrapper

## Example

```julia
# Define an ODE with Jacobian and analytical solution
function f!(du, u, p, t)
    du[1] = p[1] * u[1]
end

function jac!(J, u, p, t)
    J[1,1] = p[1]
end

function analytic(u0, p, t)
    [u0[1] * exp(p[1] * t)]
end

odefun = ODEFunction(f!; jac=jac!, analytic=analytic)
prob = ODEProblem(odefun, [1.0], (0.0, 10.0), [0.5])

# The solver can now use the Jacobian and compare with analytical solution
```

## See Also

- [Problem Interface](@ref) for using functions in problems
- SciMLBase.jl documentation for complete function specifications
- Individual solver packages for solver-specific requirements