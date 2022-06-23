# DiffEqBase.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](http://docs.sciml.ai)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai)

[![codecov](https://codecov.io/gh/SciML/DiffEqBase.jl/branch/master/graph/badge.svg?token=FwXaKBNW67)](https://codecov.io/gh/SciML/DiffEqBase.jl)
[![Build Status](https://github.com/SciML/DiffEqBase.jl/workflows/CI/badge.svg)](https://github.com/SciML/DiffEqBase.jl/actions?query=workflow%3ACI)
[![Build status](https://badge.buildkite.com/e0ee4d9d914eb44a43c291d78c53047eeff95e7edb7881b6f7.svg)](https://buildkite.com/julialang/diffeqbase-dot-jl)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)
[![Package Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/DiffEqBase)](https://pkgs.genieframework.com?packages=DiffEqBase)

DiffEqBase.jl is a component package in the [SciML Scientific Machine Learning ecosystem](https://sciml.ai/). 
It holds the sensitivity analysis utilities. Users interested in using this
functionality should check out [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl).
DiffEqBase.jl is a component package in the DiffEq ecosystem. It holds the
common types and utility functions which are shared by other component packages
in order to reduce the size of dependencies. This is so that the packages for the common interface do not require one another, allowing users to use the functionality of individual packages if they so please. Users interested in using this
functionality in full should check out [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl)

The documentation for the interfaces here can be found in [DiffEqDocs.jl](https://diffeq.sciml.ai/dev/) and [DiffEqDevDocs.jl](https://devdocs.sciml.ai/dev).
