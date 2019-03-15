# DE functions

## Function types

Functions for the different categories of differential equations are
defined by subtypes of `AbstractDiffEqFunction`:

```@docs
AbstractDiffEqFunction
AbstractDiscreteFunction
DiscreteFunction
AbstractODEFunction
ODEFunction
DynamicalODEFunction
SplitFunction
AbstractSDEFunction
SDEFunction
SplitSDEFunction
AbstractRODEFunction
RODEFunction
AbstractDDEFunction
DDEFunction
AbstractDAEFunction
DAEFunction
```

## Utilities

```@doc
has_analytic
has_jac
has_tgrad
has_invW
has_invW_t
has_paramjac
has_syms
```
