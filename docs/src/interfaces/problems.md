# Problem Interface

```@meta
CurrentModule = DiffEqBase
```

DiffEqBase extends the problem types defined in SciMLBase with additional functionality. The problem types represent the mathematical formulation of differential equations to be solved.

## Problem Types

The following abstract problem types are imported from SciMLBase and used throughout the ecosystem:

- `AbstractODEProblem`: Ordinary differential equations
- `AbstractSDEProblem`: Stochastic differential equations  
- `AbstractDDEProblem`: Delay differential equations
- `AbstractDAEProblem`: Differential-algebraic equations
- `AbstractRODEProblem`: Random ordinary differential equations
- `AbstractSDDEProblem`: Stochastic delay differential equations
- `AbstractBVProblem`: Boundary value problems
- `AbstractDiscreteProblem`: Discrete-time problems
- `AbstractNonlinearProblem`: Nonlinear equation problems
- `AbstractSteadyStateProblem`: Steady state problems
- `AbstractJumpProblem`: Jump processes
- `AbstractEnsembleProblem`: Ensemble/Monte Carlo problems

## Problem Interface

All problem types follow a common interface:

### Construction

Problems are typically constructed with:
- A function defining the dynamics
- Initial/boundary conditions
- Time span or domain
- Parameters (optional)
- Callbacks (optional)

### Common Methods

- `remake(prob; kwargs...)`: Create a new problem with modified fields
- `isinplace(prob)`: Check if the problem uses in-place operations
- `has_analytic(prob)`: Check if an analytic solution is available

## Example

```julia
# Define an ODE problem
function f(du, u, p, t)
    du[1] = p[1] * u[1]
end

u0 = [1.0]
tspan = (0.0, 10.0)
p = [0.5]

prob = ODEProblem(f, u0, tspan, p)
```

## See Also

- [Function Interface](@ref) for defining the dynamics
- [Solution Interface](@ref) for working with solutions
- SciMLBase.jl documentation for detailed problem specifications