using DiffEqBase, Test

function f(du,u,p,t)
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

@test_broken @inferred ODEProblem{true}(f,u0,tspan)
@test_broken @inferred ODEProblem(f,u0,tspan)

function f(dv,u,v,p,t)
  dv .= 2.*v
end
u0 = ones(2)
v0 = ones(2)
tspan = (0,1.0)
prob = SecondOrderODEProblem(f,u0,v0,tspan)

function f(du,u,p,t)
  du[1] = 0.2u[1]
  du[2] = 0.4u[2]
end
function g(du,u,p,t)
  du[1] = 0.2u[1]
  du[2] = 0.4u[2]
end
u0 = ones(2)
tspan = (0,1.0)
prob = SDEProblem(f,g,u0,tspan)
prob = SDEProblem{true}(f,g,u0,tspan)

@test_broken @inferred SDEProblem(f,g,u0,tspan)
@test_broken @inferred SDEProblem{true}(f,g,u0,tspan)

f_1delay = function (du,u,h,p,t)
    du[1] = - h(t-1)[1]
end
prob =  DDEProblem(f_1delay,ones(1),t->zeros(1),(0.0, 10.0),constant_lags = ones(1))
prob =  DDEProblem{true}(f_1delay,ones(1),t->zeros(1),(0.0, 10.0),dependent_lags = ones(1))

@test_broken @inferred DDEProblem(f_1delay,ones(1),t->zeros(1),(0.0, 10.0),constant_lags = ones(1))
@test_broken @inferred DDEProblem{true}(f_1delay,ones(1),t->zeros(1),(0.0, 10.0),dependent_lags = ones(1))

function f(r, yp, y, p,tres)
    r[1]  = -0.04*y[1] + 1.0e4*y[2]*y[3]
    r[2]  = -r[1] - 3.0e7*y[2]*y[2] - yp[2]
    r[1] -=  yp[1]
    r[3]  =  y[1] + y[2] + y[3] - 1.0
end
u0 = [1.0, 0, 0]
du0 = [-0.04, 0.04, 0.0]
prob_dae_resrob = DAEProblem(f,du0,u0,(0.0,100000.0))
prob_dae_resrob = DAEProblem{true}(f,du0,u0,(0.0,100000.0))

@test_broken @inferred DAEProblem(f,du0,u0,(0.0,100000.0))
@test_broken @inferred DAEProblem{true}(f,du0,u0,(0.0,100000.0))

f(u,t,W) = 1.01u.+0.87u.*W
u0 = 1.00
tspan = (0.0,1.0)
prob = RODEProblem(f,u0,tspan)
prob = RODEProblem{false}(f,u0,tspan)

@test_broken @inferred RODEProblem(f,u0,tspan)
@test_broken @inferred RODEProblem{false}(f,u0,tspan)

DiscreteProblem(ones(1),tspan)
f(t,u) = 0.5
DiscreteProblem{false}(f,ones(1),tspan)

@test_broken @inferred DiscreteProblem(ones(1),tspan)
@test_broken @inferred DiscreteProblem{false}(f,ones(1),tspan)

function f(du,u,t)
  du[1] = 2 - 2u[1]
  du[2] = u[1] - 4u[2]
end
u0 = zeros(2)
prob = SteadyStateProblem(f,u0)

@test_broken @inferred SteadyStateProblem(f,u0)
