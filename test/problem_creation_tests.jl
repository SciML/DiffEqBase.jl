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

function f(t,u,du)
  du[1] = 0.2u[1]
  du[2] = 0.4u[2]
end
function g(t,u,du)
  du[1] = 0.2u[1]
  du[2] = 0.4u[2]
end
u0 = ones(2)
tspan = (0,1.0)
prob = SDEProblem(f,g,u0,tspan)
prob = SDEProblem{true}(f,g,u0,tspan)

f_1delay = function (t,u,h,du)
    du[1] = - h(t-1)[1]
end
prob =  ConstantLagDDEProblem(f_1delay,t->zeros(1),ones(1),ones(1),(0.0, 10.0))
prob =  ConstantLagDDEProblem{true}(f_1delay,t->zeros(1),ones(1),ones(1),(0.0, 10.0))


function f(tres, y, yp, r)
    r[1]  = -0.04*y[1] + 1.0e4*y[2]*y[3]
    r[2]  = -r[1] - 3.0e7*y[2]*y[2] - yp[2]
    r[1] -=  yp[1]
    r[3]  =  y[1] + y[2] + y[3] - 1.0
end
u0 = [1.0, 0, 0]
du0 = [-0.04, 0.04, 0.0]
prob_dae_resrob = DAEProblem(f,u0,du0,(0.0,100000.0))
prob_dae_resrob = DAEProblem{true}(f,u0,du0,(0.0,100000.0))


f(t,u,W) = 1.01u.+0.87u.*W
u0 = 1.00
tspan = (0.0,1.0)
prob = RODEProblem(f,u0,tspan)
prob = RODEProblem{false}(f,u0,tspan)

DiscreteProblem(ones(1),tspan)
f(t,u) = 0.5
DiscreteProblem{false}(f,ones(1),tspan)

function f(t,u,du)
  du[1] = 2 - 2u[1]
  du[2] = u[1] - 4u[2]
end
u0 = zeros(2)
prob = SteadyStateProblem(f,u0)
