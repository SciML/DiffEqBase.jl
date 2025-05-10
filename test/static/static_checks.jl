using DiffEqBase, ComponentArrays, AllocCheck

u = ComponentArray(x=1.0, y=0.0, z=0.0)
t = 0.0
@test @allocated(DiffEqBase.ODE_DEFAULT_NORM(u, t)) < 20