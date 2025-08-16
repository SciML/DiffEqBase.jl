# Algorithm Interface

```@meta
CurrentModule = DiffEqBase
```

Algorithms define the numerical methods used to solve differential equation problems. DiffEqBase provides the interface that all algorithm implementations must follow.

## Algorithm Types

The following abstract algorithm types are imported from SciMLBase:

- `AbstractDEAlgorithm`: Base type for all DE algorithms
- `AbstractODEAlgorithm`: Algorithms for ODEs
- `AbstractSDEAlgorithm`: Algorithms for SDEs
- `AbstractDDEAlgorithm`: Algorithms for DDEs
- `AbstractDAEAlgorithm`: Algorithms for DAEs
- `AbstractRODEAlgorithm`: Algorithms for RODEs
- `AbstractSDDEAlgorithm`: Algorithms for SDDEs
- `AbstractBVPAlgorithm`: Algorithms for BVPs
- `AbstractSteadyStateAlgorithm`: Algorithms for steady states
- `AbstractNonlinearAlgorithm`: Algorithms for nonlinear problems
- `DAEInitializationAlgorithm`: Algorithms for DAE initialization
- `AbstractSensitivityAlgorithm`: Algorithms for sensitivity analysis
- `EnsembleAlgorithm`: Algorithms for ensemble problems

## Algorithm Interface

### Properties

Algorithms can have various properties that can be queried:

- `isadaptive(alg)`: Whether the algorithm uses adaptive time-stepping
- `isdiscrete(alg)`: Whether the algorithm is for discrete problems
- `isautodifferentiable(alg)`: Whether the algorithm supports AD

### Common Algorithm Options

When solving with an algorithm, common options include:

- `abstol`: Absolute tolerance
- `reltol`: Relative tolerance  
- `dtmax`: Maximum time step
- `dtmin`: Minimum time step
- `saveat`: Times to save the solution
- `save_everystep`: Whether to save all steps
- `callback`: Callbacks to apply
- `maxiters`: Maximum iterations

## Ensemble Algorithms

Special algorithms for ensemble/Monte Carlo simulations:

- `EnsembleThreads()`: Thread-based parallelism
- `EnsembleDistributed()`: Distributed parallelism
- `EnsembleSerial()`: Serial execution

## Creating Custom Algorithms

To create a custom algorithm:

```julia
struct MyAlgorithm <: AbstractODEAlgorithm
    # Algorithm parameters
end

# Implement required methods
isadaptive(::MyAlgorithm) = false
```

## Example

```julia
# Use a specific algorithm
using OrdinaryDiffEq

prob = ODEProblem((u,p,t) -> -u, 1.0, (0.0, 10.0))

# Solve with Tsit5 algorithm
sol = solve(prob, Tsit5(); abstol=1e-8, reltol=1e-8)

# Check algorithm properties
isadaptive(Tsit5())  # true
```

## See Also

- [Problem Interface](@ref) for problem types
- [Solution Interface](@ref) for solution types
- Individual solver packages for specific algorithms
- SciMLBase.jl documentation for complete algorithm specifications