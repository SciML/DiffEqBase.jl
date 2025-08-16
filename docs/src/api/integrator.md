# Integrator Interface

```@meta
CurrentModule = DiffEqBase
```

The integrator interface provides low-level control over the integration process. This allows for step-by-step integration and fine-grained control over the solution process.

## Integrator Accessors

These functions provide access to and modification of integrator state:

```@autodocs
Modules = [DiffEqBase]
Pages = ["integrator_accessors.jl"]
```

## Internal Integration Methods

### Euler Method

The internal Euler method implementation for basic integration:

```@autodocs
Modules = [DiffEqBase]
Pages = ["internal_euler.jl"]
```

### Interpolation

Internal interpolation utilities:

```@autodocs
Modules = [DiffEqBase]
Pages = ["internal_itp.jl"]
```