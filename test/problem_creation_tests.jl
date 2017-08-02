using DiffEqBase, Base.Test

function f(t,u,du)
  du[1] = 0.2u[1]
  du[2] = 0.4u[2]
end
u0 = ones(2)
tspan = (0,1.0)

prob = ODEProblem(f,u0,tspan)
@test typeof(prob.tspan) == Tuple{Float64,Float64}
prob = ODEProblem{true}(f,u0,tspan)
@test typeof(prob.tspan) == Tuple{Float64,Float64}
@test isinplace(prob) == true
prob = ODEProblem{false}(f,u0,tspan)
@test isinplace(prob) == false

function f(t,u,v,dv)
  dv .= 2.*v
end
u0 = ones(2)
v0 = ones(2)
tspan = (0,1.0)
prob = SecondOrderODEProblem(f,u0,v0,tspan)
