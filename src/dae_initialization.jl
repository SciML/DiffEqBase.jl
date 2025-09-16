# DAE Initialization Algorithms
# The base types are in SciMLBase, we just add the specific algorithms here

import SciMLBase: DAEInitializationAlgorithm, initialize_dae!, CheckInit, NoInit, OverrideInit

# Re-export the SciMLBase initialization algorithms for convenience
# Note: Docstrings for these types should be added in SciMLBase where they are defined
export CheckInit, NoInit, OverrideInit

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