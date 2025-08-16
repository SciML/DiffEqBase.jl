# Solution Interface

```@meta
CurrentModule = DiffEqBase
```

Solutions represent the results of solving differential equation problems. DiffEqBase extends the solution types from SciMLBase with additional functionality.

## Solution Types

The following abstract solution types are imported from SciMLBase:

- `AbstractODESolution`: Solutions to ODE problems
- `AbstractRODESolution`: Solutions to RODE problems
- `AbstractDAESolution`: Solutions to DAE problems
- `AbstractDDESolution`: Solutions to DDE problems
- `AbstractEnsembleSolution`: Solutions to ensemble problems
- `NonlinearSolution`: Solutions to nonlinear problems
- `AbstractTimeseriesSolution`: Base type for time-dependent solutions
- `AbstractNoTimeSolution`: Base type for non-time-dependent solutions
- `AbstractAnalyticalSolution`: Solutions with analytical expressions

## Solution Interface

### Accessing Solution Data

Solutions typically contain:
- `sol.t`: Time points
- `sol.u`: Solution values at time points
- `sol.prob`: The original problem
- `sol.alg`: The algorithm used
- `sol.retcode`: Return code indicating success/failure
- `sol.stats`: Statistics about the solution process

### Interpolation

Solutions support interpolation between computed points:

```julia
sol(t)  # Interpolate solution at time t
```

### Common Methods

- `sol[i]`: Access i-th solution point
- `length(sol)`: Number of solution points
- `size(sol)`: Dimensions of the solution
- `plot(sol)`: Visualize the solution (requires Plots.jl)

## Statistics

Solution statistics can include:
- Number of function evaluations
- Number of Jacobian evaluations
- Number of accepted/rejected steps
- Step size information

Access via `sol.stats`.

## Error Analysis

For problems with known analytical solutions:

```julia
calculate_solution_errors!(sol)
```

This computes various error metrics when an analytical solution is available.

## Ensemble Solutions

Ensemble solutions contain multiple trajectories:

```julia
ensemble_sol[i]  # Access i-th trajectory
calculate_ensemble_errors(ensemble_sol)  # Compute ensemble statistics
```

## Example

```julia
# Solve an ODE and work with the solution
prob = ODEProblem((u,p,t) -> -u, 1.0, (0.0, 10.0))
sol = solve(prob)

# Access solution
sol.t        # Time points
sol.u        # Solution values
sol(5.0)     # Interpolate at t=5.0
sol.stats    # Solution statistics
```

## See Also

- [Problem Interface](@ref) for creating problems
- [Statistics](@ref) for detailed statistics functions
- SciMLBase.jl documentation for complete solution specifications