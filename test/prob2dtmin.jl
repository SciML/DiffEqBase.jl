using DiffEqBase, ForwardDiff, Test
@test DiffEqBase.prob2dtmin((0.0, 10.0), 1.0, false) == eps(Float64)
@test DiffEqBase.prob2dtmin((0f0, 10f0), 1f0, false) == eps(Float32)
@test DiffEqBase.prob2dtmin((0.0, 10.0), ForwardDiff.Dual(1.0), false) == eps(Float64)
