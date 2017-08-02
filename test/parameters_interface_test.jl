using OrdinaryDiffEq
using ParameterizedFunctions
using DiffEqBase
using Base.Test

# sweet sweet lorenz equation
lorenz = @ode_def Lorenz begin
  dx = σ*(y-x)
  dy = ρ*x-y-x*z
  dz = x*y-β*z
end σ => 10.0 β = 8.0/3.0 ρ => 28.0

lorenz_compat = DiffEqBase.CompatibilityDiffEqFunction(lorenz)
@test lorenz_compat.jac != nothing

# make a clean instance to change
lorenz_test = Lorenz()
lorenz_test.params

# make sure we need to provide all the free parameters
@test_throws DimensionMismatch set_param_values!(lorenz_test, [0.5, 3.2, 5.8])
@test_throws DimensionMismatch set_param_values!(lorenz_test, [0.5])

# make sure the general setter works
set_param_values!(lorenz_test, [0.5, 3.2])
@test lorenz_test.σ == 0.5
@test lorenz_test.ρ == 3.2

set_param_values!(lorenz_test, Dict(:ρ => 11.3))
@test lorenz_test.ρ == 11.3
@test lorenz_test.σ == 0.5
set_param_values!(lorenz_test, Dict(:σ => -2.1, :ρ => 5.3211115))
@test lorenz_test.σ == -2.1
@test lorenz_test.ρ == 5.3211115

# This should throw an error, not sure how to check for it
#set_param_values!(lorenz_test, Dict(:σ => 6.0, :ρ => 4.3, :β => 1.1))
