using Plots, OrdinaryDiffEq
unicodeplots()

f(x,p,t) = p*x
prob = ODEProblem(f,[1.0,2.0,3.0],(0.0,1.0),-.2)
sol = solve(prob, Tsit5())

f1(t,x,y) = (t,x+y)
f2(t,x) = (t,x)
f3(t,x,y) = (t,x)
plot(sol, vars = [(f1,0,1,2),])
plot(sol, vars = [(f1,0,1,2), (0,3)])
plot(sol, vars = [(f1,0,1,2), (f2,0,3)])
plot(sol, vars = [(f1,0,1,2), (f3,0,3,1)])
