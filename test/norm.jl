using Test
using ForwardDiff: Dual

using DiffEqBase: ODE_DEFAULT_NORM
const internalnorm = ODE_DEFAULT_NORM

val = rand(10)
par = rand(10)
u = Dual.(val, par)
@test internalnorm(val, 1) == internalnorm(u, 1)
