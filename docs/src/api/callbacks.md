# Callbacks

```@meta
CurrentModule = DiffEqBase
```

Callbacks provide a way to affect the integrator during the solution process. They can be used to implement event handling, saving, and control mechanisms.

## Initialization and Finalization

```@docs
initialize!
finalize!
```

## Callback Utilities

The callback system integrates with the SciMLBase callback infrastructure. The main types are:

- `CallbackSet`: A collection of callbacks
- `ContinuousCallback`: Callbacks triggered by continuous conditions
- `DiscreteCallback`: Callbacks triggered at discrete time points

These base types are defined in SciMLBase.jl and extended here with initialization and finalization methods.

## Helper Functions

```@autodocs
Modules = [DiffEqBase]
Pages = ["callbacks.jl"]
Filter = t -> !occursin("initialize!", string(t)) && !occursin("finalize!", string(t))
```