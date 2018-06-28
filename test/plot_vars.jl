using DiffEqBase.InternalEuler, DiffEqBase, Test

# Here's the problem to solve

struct LorenzFunction <: Function
    syms::Vector{Symbol}
end

function (::LorenzFunction)(u,p,t)
 [10.0(u[2]-u[1]),u[1]*(28.0-u[3]) - u[2],u[1]*u[2] - (8/3)*u[3]]
end
lorenz = LorenzFunction([:x,:y,:z])

u0 = [1., 5., 10.]
tspan = (0., 100.)
prob = ODEProblem(lorenz, u0, tspan)
dt = 0.1
sol = solve(prob,InternalEuler.FwdEulerAlg(), tstops=0:dt:1)

@test DiffEqBase.interpret_vars([(0,1), (1,3), (4,5)],sol) == [(DiffEqBase.DEFAULT_PLOT_FUNC,0,1), (DiffEqBase.DEFAULT_PLOT_FUNC,1,3), (DiffEqBase.DEFAULT_PLOT_FUNC,4,5)]
@test DiffEqBase.interpret_vars([1, (1,3), (4,5)],sol) == [(DiffEqBase.DEFAULT_PLOT_FUNC,0,1), (DiffEqBase.DEFAULT_PLOT_FUNC,1,3), (DiffEqBase.DEFAULT_PLOT_FUNC,4,5)]
@test DiffEqBase.interpret_vars([1, 3, 4],sol) == [(DiffEqBase.DEFAULT_PLOT_FUNC,0,1), (DiffEqBase.DEFAULT_PLOT_FUNC,0,3), (DiffEqBase.DEFAULT_PLOT_FUNC,0,4)]
@test DiffEqBase.interpret_vars(([1,2,3], [4,5,6]),sol) == [(DiffEqBase.DEFAULT_PLOT_FUNC,1,4), (DiffEqBase.DEFAULT_PLOT_FUNC,2,5), (DiffEqBase.DEFAULT_PLOT_FUNC,3,6)]
@test DiffEqBase.interpret_vars((1, [2,3,4]),sol) == [(DiffEqBase.DEFAULT_PLOT_FUNC,1,2), (DiffEqBase.DEFAULT_PLOT_FUNC,1,3), (DiffEqBase.DEFAULT_PLOT_FUNC,1,4)]

@test DiffEqBase.interpret_vars([(:t,:x),(:t,:y)],sol) == [(DiffEqBase.DEFAULT_PLOT_FUNC,0,1), (DiffEqBase.DEFAULT_PLOT_FUNC,0,2)]
@test DiffEqBase.interpret_vars([:x, (0,:x), (:x,:y)],sol) == [(DiffEqBase.DEFAULT_PLOT_FUNC,0,1), (DiffEqBase.DEFAULT_PLOT_FUNC,0,1), (DiffEqBase.DEFAULT_PLOT_FUNC,1,2)]
@test DiffEqBase.interpret_vars([:x, :y, :z],sol) == [(DiffEqBase.DEFAULT_PLOT_FUNC,0,1), (DiffEqBase.DEFAULT_PLOT_FUNC,0,2), (DiffEqBase.DEFAULT_PLOT_FUNC,0,3)]
@test DiffEqBase.interpret_vars(([:x,:x], [:y,:z]),sol) == [(DiffEqBase.DEFAULT_PLOT_FUNC,1,2), (DiffEqBase.DEFAULT_PLOT_FUNC,1,3)]
@test DiffEqBase.interpret_vars((:x, [:y,:z]),sol) == [(DiffEqBase.DEFAULT_PLOT_FUNC,1,2), (DiffEqBase.DEFAULT_PLOT_FUNC,1,3)]

f(x,y) = (x+y,y)
@test DiffEqBase.interpret_vars([(f,0,1), (1,3), (4,5)],sol) == [(f,0,1), (DiffEqBase.DEFAULT_PLOT_FUNC,1,3), (DiffEqBase.DEFAULT_PLOT_FUNC,4,5)]
@test DiffEqBase.interpret_vars([1, (f,1,3), (4,5)],sol) == [(DiffEqBase.DEFAULT_PLOT_FUNC,0,1), (f,1,3), (DiffEqBase.DEFAULT_PLOT_FUNC,4,5)]
@test DiffEqBase.interpret_vars([(f,:t,:x),(:t,:y)],sol) == [(f,0,1), (DiffEqBase.DEFAULT_PLOT_FUNC,0,2)]
@test DiffEqBase.interpret_vars([:x, (f,0,:x), (:x,:y)],sol) == [(DiffEqBase.DEFAULT_PLOT_FUNC,0,1), (f,0,1), (DiffEqBase.DEFAULT_PLOT_FUNC,1,2)]
@test DiffEqBase.interpret_vars([(:x,:y)],sol) == [(DiffEqBase.DEFAULT_PLOT_FUNC,1,2)]
@test DiffEqBase.interpret_vars((f,:x,:y),sol) == [(f,1,2)]
