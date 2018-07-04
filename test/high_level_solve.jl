using DiffEqBase, Test

@test DiffEqBase.promote_tspan((0.0,1.0)) == (0.0,1.0)
@test DiffEqBase.promote_tspan((0,1.0)) == (0.0,1.0)
@test DiffEqBase.promote_tspan(1.0) == (0.0,1.0)
@test DiffEqBase.promote_tspan(nothing) == (nothing,nothing)
