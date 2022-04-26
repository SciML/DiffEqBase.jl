using OrdinaryDiffEq, Test

f(u,p,t) = 2u
u0 = 0.5
tspan = (0.0,1.0)
prob = ODEProblem(f,g,u0,tspan)
sol = solve(prob,Tsit5())

function f(du, u, p, t)
    du .= 2.0 * u
end
sol = solve(prob,Tsit5())