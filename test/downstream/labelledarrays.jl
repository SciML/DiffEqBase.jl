using OrdinaryDiffEq
using LabelledArrays

function f(out,du,u,p,t)
  out.x = - 0.04u.x              + 1e4*u.y*u.z - du.x
  out.y = + 0.04u.x - 3e7*u.y^2 - 1e4*u.y*u.z - du.y
  out.z = u.x + u.y + u.z - 1.0
end

u₀ = LVector(x=1.0, y=0.0, z=0.0)
du₀ = LVector(x=-0.04, y=0.04, z=0.0)
tspan = (0.0,100000.0)

differential_vars = LVector(x=true, y=true, z=false)
prob = DAEProblem(f,du₀,u₀,tspan,differential_vars=differential_vars)

sol = solve(prob, DImplicitEuler())
