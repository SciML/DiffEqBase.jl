# Problems


TODO:
* How do optional parameters (`p=nothing`) work?


## Problem types

```@docs
DEProblem
```

### Discrete problems

```@docs
AbstractDiscreteProblem
DiscreteProblem
```

### ODE problems

```@docs
AbstractODEProblem
ODEProblem
```

#### Subtypes

TODO: document `problem_type` field.

```@docs
StandardODEProblem
AbstractDynamicalODEProblem
DynamicalODEProblem
SecondOrderODEProblem
AbstractSplitODEProblem
SplitODEProblem
```

### Steady state problems

```@docs
AbstractSteadyStateProblem
SteadyStateProblem
```

### SDE problems

```@docs
AbstractSDEProblem
```

### RODE problems

```@docs
AbstractRODEProblem
```

### DDE problems

```@docs
AbstractDDEProblem
```

### DAE problems

```@docs
AbstractDAEProblem
```

### Jump problems

```@docs
AbstractJumpProblem
```
