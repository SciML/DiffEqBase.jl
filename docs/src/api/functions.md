# DE functions

## Function types

Functions for the different categories of differential equations are
defined by subtypes of `AbstractDiffEqFunction`:

```@docs
DiffEqBase.AbstractDiffEqFunction
DiffEqBase.AbstractDiscreteFunction
DiscreteFunction
DiffEqBase.AbstractODEFunction
ODEFunction
DynamicalODEFunction
SplitFunction
DiffEqBase.AbstractSDEFunction
SDEFunction
SplitSDEFunction
DiffEqBase.AbstractRODEFunction
RODEFunction
DiffEqBase.AbstractDDEFunction
DDEFunction
DiffEqBase.AbstractDAEFunction
DAEFunction
```

## Utilities

```@docs
has_analytic
has_jac
has_tgrad
has_invW
has_invW_t
has_paramjac
has_syms
```
