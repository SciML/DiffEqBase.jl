# DiffEqBase.jl

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://github.com/SciML/DiffEqBase.jl/workflows/CI/badge.svg)](https://github.com/SciML/DiffEqBase.jl/actions?query=workflow%3ACI)
[![GitlabCI](https://gitlab.com/JuliaGPU/DiffEqBase.jl/badges/master/pipeline.svg)](https://gitlab.com/JuliaGPU/DiffEqBase.jl/pipelines)

DiffEqBase.jl is a component package in the DiffEq ecosystem. It holds the
common types and utility functions which are shared by other component packages
in order to reduce the size of dependencies. This is so that the packages for the common interface do not require one another, allowing users to use the functionality of individual packages if they so please. Users interested in using this
functionality in full should check out [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl)

The documentation for the interfaces here can be found in [DiffEqDocs.jl](https://diffeq.sciml.ai/dev/) and [DiffEqDevDocs.jl](https://devdocs.sciml.ai/dev). Specific parts to note are:
