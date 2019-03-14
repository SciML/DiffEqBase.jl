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

### BVP problems

```@docs
AbstractBVPProblem
StandardBVProblem
BVProblem
TwoPointBVProblem
```

### SDE problems

```@docs
AbstractSDEProblem
StandardSDEProblem
SDEProblem
AbstractSplitSDEProblem
SplitSDEProblem
```

### RODE problems

```@docs
AbstractRODEProblem
RODEProblem
```

### DDE problems

```@docs
AbstractDDEProblem
DDEProblem
```

### DAE problems

```@docs
AbstractDAEProblem
DAEProblem
```

### Jump problems

```@docs
AbstractJumpProblem
```
