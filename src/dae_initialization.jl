# DAE Initialization Algorithms
# The base types are in SciMLBase, we just add the specific algorithms here

import SciMLBase: DAEInitializationAlgorithm, initialize_dae!, CheckInit, NoInit, OverrideInit

# Re-export the SciMLBase initialization algorithms for convenience
export CheckInit, NoInit, OverrideInit

"""
    struct CheckInit <: DAEInitializationAlgorithm

An initialization algorithm that only checks if the initial conditions are consistent
with the DAE constraints, without attempting to modify them. If the conditions are not
consistent within the solver's tolerance, an error will be thrown.

This is useful when:
- You have already computed consistent initial conditions
- You want to verify the consistency of your initial guess
- You want to ensure no automatic modifications are made to your initial conditions

## Example
```julia
prob = DAEProblem(f, du0, u0, tspan)
sol = solve(prob, IDA(), initializealg = CheckInit())
```
"""
CheckInit

"""
    struct NoInit <: DAEInitializationAlgorithm

An initialization algorithm that completely skips the initialization phase. The solver
will use the provided initial conditions directly without any consistency checks or
modifications.

⚠️ **Warning**: Using `NoInit()` with inconsistent initial conditions will likely cause
solver failures or incorrect results. Only use this when you are absolutely certain
your initial conditions satisfy all DAE constraints.

This is useful when:
- You know your initial conditions are already perfectly consistent
- You want to avoid the computational cost of initialization
- You are debugging solver issues and want to isolate initialization from integration

## Example
```julia
prob = DAEProblem(f, du0_consistent, u0_consistent, tspan)
sol = solve(prob, IDA(), initializealg = NoInit())
```
"""
NoInit

"""
    struct OverrideInit <: DAEInitializationAlgorithm

An initialization algorithm that uses a separate initialization problem to find
consistent initial conditions. This is typically used with ModelingToolkit.jl
which can generate specialized initialization problems based on the model structure.

When using `OverrideInit`, the problem must have `initialization_data` that contains
an `initializeprob` field with the initialization problem to solve.

This algorithm is particularly useful for:
- High-index DAEs that have been index-reduced
- Systems with complex initialization requirements
- ModelingToolkit models with custom initialization equations

## Example
```julia
# Typically used automatically with ModelingToolkit
@named sys = ODESystem(eqs, t, vars, params)
sys = structural_simplify(dae_index_lowering(sys))
prob = DAEProblem(sys, [], (0.0, 1.0), [])
# Will automatically use OverrideInit if initialization_data exists
sol = solve(prob, IDA())
```
"""
OverrideInit

"""
    struct DefaultInit <: DAEInitializationAlgorithm

The default initialization algorithm for DAEs. This will use heuristics to
determine the most appropriate initialization based on the problem type.

For Sundials, this will use:
- `OverrideInit` if the problem has `initialization_data` (typically from ModelingToolkit)
- `CheckInit` otherwise
"""
struct DefaultInit <: DAEInitializationAlgorithm end

"""
    struct BrownBasicInit <: DAEInitializationAlgorithm

The Brown basic initialization algorithm for DAEs. This implementation
is based on the algorithm described in:

Peter N. Brown, Alan C. Hindmarsh, and Linda R. Petzold,
"Consistent Initial Condition Calculation for Differential-Algebraic Systems",
SIAM Journal on Scientific Computing, Vol. 19, No. 5, pp. 1495-1512, 1998.
DOI: https://doi.org/10.1137/S1064827595289996

This method modifies the algebraic variables and their derivatives to be
consistent with the DAE constraints, while keeping the differential variables
fixed. It uses Newton's method to solve for consistent initial values.

This is the default initialization for many DAE solvers when `differential_vars`
is provided, allowing the solver to distinguish between differential and algebraic
variables.
"""
struct BrownBasicInit <: DAEInitializationAlgorithm end

"""
    struct ShampineCollocationInit <: DAEInitializationAlgorithm

The Shampine collocation initialization algorithm for DAEs. This implementation
is based on the algorithm described in:

Lawrence F. Shampine, "Consistent Initial Condition for Differential-Algebraic
Systems", SIAM Journal on Scientific Computing, Vol. 22, No. 6, pp. 2007-2026, 2001.
DOI: https://doi.org/10.1137/S1064827599355049

This method uses collocation on the first two steps to find consistent initial
conditions. It modifies both the differential and algebraic variables to satisfy
the DAE constraints. This is more general than BrownBasicInit but may be more
expensive computationally.

This method is useful when you need to modify all variables (both differential
and algebraic) to achieve consistency, rather than just the algebraic ones.
"""
struct ShampineCollocationInit <: DAEInitializationAlgorithm end

export DefaultInit, BrownBasicInit, ShampineCollocationInit