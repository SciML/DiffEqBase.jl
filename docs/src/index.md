# DiffEqBase.jl

```@meta
CurrentModule = DiffEqBase
```

## Overview

DiffEqBase.jl is the common foundation library for the DifferentialEquations.jl ecosystem. It provides the core types, interfaces, and functionality that are shared across all differential equation solvers in the SciML ecosystem.

## Features

- Common problem and solution types for ODEs, SDEs, DDEs, DAEs, and more
- Unified integrator interface for all differential equation types
- Callback system for event handling and control
- Statistics and analysis utilities
- Interpolation and tableau definitions
- Integration with automatic differentiation and sensitivity analysis

## Installation

```julia
using Pkg
Pkg.add("DiffEqBase")
```

## Getting Started

DiffEqBase.jl is typically used as a dependency by other packages in the DifferentialEquations.jl ecosystem. However, you can use it directly for low-level operations:

```julia
using DiffEqBase

# Most functionality is re-exported from SciMLBase
# and used through higher-level packages
```

## Main Components

### Problems and Solutions
The package defines abstract types and interfaces for various differential equation problems and their solutions.

### Integrators
Low-level integrator interface for step-by-step integration control.

### Callbacks
Event handling system for modifying solutions during integration.

### Utilities
Various utility functions for working with differential equations.

## Package Ecosystem

DiffEqBase.jl is part of the larger SciML ecosystem:

- **SciMLBase.jl**: Core types and interfaces
- **OrdinaryDiffEq.jl**: ODE solvers
- **StochasticDiffEq.jl**: SDE solvers
- **DelayDiffEq.jl**: DDE solvers
- **DifferentialEquations.jl**: Meta-package

## Index

```@index
```