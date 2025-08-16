# Internal Functions

```@meta
CurrentModule = DiffEqBase
```

Internal functions that are part of the public API but primarily used by solver packages.

## Solver Functions

```@autodocs
Modules = [DiffEqBase]
Pages = ["solve.jl"]
```

## No-Recompile Utilities

Functions for avoiding recompilation:

```@autodocs
Modules = [DiffEqBase]
Pages = ["norecompile.jl"]
```

## Cost Functions

```@docs
DECostFunction
```

## Parameterized Functions

```@docs
AbstractParameterizedFunction
```

## Convergence Setup

```@docs
ConvergenceSetup
```

## Sensitivity Analysis

```@docs
SensitivityADPassThrough
```

## Keyword Argument Handling

```@docs
KeywordArgError
KeywordArgWarn
KeywordArgSilent
```