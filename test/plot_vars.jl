using OrdinaryDiffEq, ParameterizedFunctions, DiffEqBase, Base.Test

# Here's the problem to solve

lorenz = @ode_def Lorenz begin
  dx = σ*(y-x)
  dy = ρ*x-y-x*z
  dz = x*y-β*z
end σ = 10. β = 8./3. ρ => 28.

u0 = [1., 5., 10.]
tspan = (0., 100.)
prob = ODEProblem(lorenz, u0, tspan)
sol = solve(prob,Tsit5())

@test DiffEqBase.interpret_vars([(0,1), (1,3), (4,5)],sol) == [(0,1), (1,3), (4,5)]
@test DiffEqBase.interpret_vars([1, (1,3), (4,5)],sol) == [(0,1), (1,3), (4,5)]
@test DiffEqBase.interpret_vars([1, 3, 4],sol) == [(0,1), (0,3), (0,4)]
@test DiffEqBase.interpret_vars(([1,2,3], [4,5,6]),sol) == [(1,4), (2,5), (3,6)]
@test DiffEqBase.interpret_vars((1, [2,3,4]),sol) == [(1,2), (1,3), (1,4)]

@test DiffEqBase.interpret_vars([(:t,:x),(:t,:y)],sol) == [(0,1), (0,2)]
@test DiffEqBase.interpret_vars([:x, (0,:x), (:x,:y)],sol) == [(0,1), (0,1), (1,2)]
@test DiffEqBase.interpret_vars([:x, :y, :z],sol) == [(0,1), (0,2), (0,3)]
@test DiffEqBase.interpret_vars(([:x,:x], [:y,:z]),sol) == [(1,2), (1,3)]
@test DiffEqBase.interpret_vars((:x, [:y,:z]),sol) == [(1,2), (1,3)]
