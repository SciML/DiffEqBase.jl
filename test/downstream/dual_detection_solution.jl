using OrdinaryDiffEq

## https://github.com/SciML/DifferentialEquations.jl/issues/1013

mutable struct SomeObject
    position
    velocity
    trajectory
end

object = SomeObject(0, 1, nothing)

# Current dynamics don't involve the object for the sake of MWE, but they could. 
function dynamics(du, u, p, t)
    du[1] = u[2]
    du[2] = -u[2]
end

for i in 1:2
    initial_state = [0, 0]
    tspan = (0.0, 5.0)
    prob = ODEProblem(dynamics, initial_state, tspan, object)
    sol = solve(prob, Tsit5())
    object.trajectory = sol
end