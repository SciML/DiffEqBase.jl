using LinearAlgebra, OrdinaryDiffEq, Test
import ForwardDiff

# setup
pd = 3

## ode with complex numbers
H0 = rand(ComplexF64,pd,pd)
A = rand(ComplexF64,pd,pd)
function f!(du,u,p,t)
  a,b,c = p
  du .= (A*u) * (a*cos(b*t+c))
  du.+= H0*u
  return nothing
end

## time span
tspan = (0.0, 1.0)

## initial state
u0 = hcat(normalize(rand(ComplexF64,pd)), normalize(rand(pd)))

## ode problem
prob0 = ODEProblem(f!, u0, tspan, rand(3); saveat=range(tspan..., 3), reltol=1e-6, alg=Tsit5())
## final state cost
cost(u) = abs2(tr(first(u)'u[2])) - abs2(tr(first(u)'last(u)))

## real loss function via complex ode
function loss(p)
  prob = remake(prob0; p)
  sol = solve(prob)
  cost(sol.u) + sum(p) / 10
end

## same problem via reals
### realify complex ode problem
function real_f(du,u,p,t)
  complex_u = complex.(selectdim(u,3,1), selectdim(u,3,2))
  complex_du = copy(complex_u)
  prob0.f(complex_du, complex_u, p, t)
  selectdim(du,3,1) .= real(complex_du)
  selectdim(du,3,2) .= imag(complex_du)
  return nothing
end
prob0_real = remake(prob0; f=real_f, u0=cat(real(prob0.u0), imag(prob0.u0); dims=3))
### real loss function via real ode
function loss_via_real(p)
  prob = remake(prob0_real; p)
  sol = solve(prob)
  u = [complex.(selectdim(u,3,1), selectdim(u,3,2)) for u=sol.u]
  cost(u) + sum(p) / 10
end

# assert
@assert eltype(last(solve(prob0     ).u)) <: Complex
@assert eltype(last(solve(prob0_real).u)) <: Real
function assert_fun()
  p0 = rand(3)
  isapprox(loss(p0), loss_via_real(p0); rtol=1e-4)
end
@assert all([assert_fun() for _=1:2^6])

# test ad with ForwardDiff
function test_ad()
  p0 = rand(3)
  grad_real = ForwardDiff.gradient(loss_via_real, p0)
  grad_complex = ForwardDiff.gradient(loss, p0)
  any(isnan.(grad_complex)) && @warn "NaN detected in gradient using ode with complex numbers !!"
  any(isnan.(grad_real   )) && @warn "NaN detected in gradient using realified ode !!"
  rel_err = norm(grad_complex-grad_real)/max(norm(grad_complex), norm(grad_real))
  isapprox(grad_complex, grad_real; rtol=1e-6) ? true : (@show rel_err; false)
end

@time @test all([test_ad() for _=1:2^6])
