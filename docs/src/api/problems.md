# Problems


TODO:
* How do optional parameters (`p=nothing`) work?

## Problem types

```@docs
DiffEqBase.DEProblem
```

### Discrete problems

```@docs
DiffEqBase.AbstractDiscreteProblem
DiscreteProblem
```

### ODE problems

```@docs
DiffEqBase.AbstractODEProblem
ODEProblem
```

#### Subtypes

TODO: document `problem_type` field.

```@docs
DiffEqBase.StandardODEProblem
DiffEqBase.AbstractDynamicalODEProblem
DynamicalODEProblem
SecondOrderODEProblem
DiffEqBase.AbstractSplitODEProblem
SplitODEProblem
```

### Steady state problems

```@docs
DiffEqBase.AbstractSteadyStateProblem
SteadyStateProblem
```

### Boundary value problems

```@docs
DiffEqBase.AbstractBVProblem
DiffEqBase.StandardBVProblem
BVProblem
TwoPointBVProblem
```

### SDE problems

```@docs
DiffEqBase.AbstractSDEProblem
DiffEqBase.StandardSDEProblem
SDEProblem
DiffEqBase.AbstractSplitSDEProblem
SplitSDEProblem
```

### RODE problems

```@docs
DiffEqBase.AbstractRODEProblem
RODEProblem
```

### DDE problems

```@docs
DiffEqBase.AbstractDDEProblem
DDEProblem
```

### DAE problems

```@docs
DiffEqBase.AbstractDAEProblem
DAEProblem
```

### Jump problems

```@docs
DiffEqBase.AbstractJumpProblem
```


## Utilities

```@docs
isinplace
DiffEqBase.is_diagonal_noise
DiffEqBase.promote_tspan
```
