# Tableaus

```@meta
CurrentModule = DiffEqBase
```

Butcher tableaus and related structures for Runge-Kutta methods.

## Tableau Types

```@docs
Tableau
ODERKTableau
```

## Tableau Functions

```@autodocs
Modules = [DiffEqBase]
Pages = ["tableaus.jl"]
Filter = t -> !(t in [Tableau, ODERKTableau])
```