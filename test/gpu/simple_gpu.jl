using OrdinaryDiffEq, CUDA, LinearAlgebra, Test
Base.:+(A::CuArray,I::UniformScaling) = A + CuArray(I,size(A,1),size(A,2))
Base.:-(A::CuArray,I::UniformScaling) = A + CuArray(I,size(A,1),size(A,2))
function f(u,p,t)
    A*u
end
function f(du,u,p,t)
    mul!(du,A,u)
end
function jac(J,u,p,t)
    J .= A
end
function jac(u,p,t)
    A
end
function tgrad(du,u,p,t)
    du .= 0
end
function tgrad(u,p,t)
    zero(u)
end
ff = ODEFunction(f,jac=jac,tgrad=tgrad)
CUDA.allowscalar(false)
A = cu(-rand(3,3))
u0 = cu([1.0;0.0;0.0])
tspan = (0.,100.)

prob = ODEProblem(ff,u0,tspan)
sol = solve(prob,Tsit5())
@test solve(prob,Rosenbrock23()).retcode == :Success
solve(prob,Rosenbrock23(autodiff=false))

prob_oop = ODEProblem{false}(ff,u0,tspan)
CUDA.allowscalar(false)
sol = solve(prob_oop,Tsit5())
@test solve(prob_oop,Rosenbrock23()).retcode == :Success
@test solve(prob_oop,Rosenbrock23(autodiff=false)).retcode == :Success

prob_nojac = ODEProblem(f,u0,tspan)
@test solve(prob_nojac,Rosenbrock23()).retcode == :Success
@test solve(prob_nojac,Rosenbrock23(autodiff=false)).retcode == :Success
@test solve(prob_nojac,Rosenbrock23(autodiff=false,diff_type = Val{:central})).retcode == :Success
@test solve(prob_nojac,Rosenbrock23(autodiff=false,diff_type = Val{:complex})).retcode == :Success

#=
prob_nojac_oop = ODEProblem{false}(f,u0,tspan)
DiffEqBase.prob2dtmin(prob_nojac_oop)
@test_broken solve(prob_nojac_oop,Rosenbrock23()).retcode == :Success
@test_broken solve(prob_nojac_oop,Rosenbrock23(autodiff=false)).retcode == :Success
@test_broken solve(prob_nojac_oop,Rosenbrock23(autodiff=false,diff_type = Val{:central})).retcode == :Success
@test_broken solve(prob_nojac_oop,Rosenbrock23(autodiff=false,diff_type = Val{:complex})).retcode == :Success
=#

# Test auto-offload
_A = -rand(3,3)
function f2(du,u,p,t)
    mul!(du,_A,u)
end
function jac2(J,u,p,t)
    J .= _A
end
ff2 = ODEFunction(f2,jac=jac2)
u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob_num = ODEProblem(ff2,u0,tspan)
sol = solve(prob_num,Rosenbrock23(linsolve=LinSolveGPUFactorize()))

# Complex Numbers Adaptivity DifferentialEquations.jl#460
f_complex(u,nothing,t) = 1/2 .*u
u0 = cu(rand(32,32).+ 1im*rand(32,32));
prob = ODEProblem(f_complex,u0,(0.0,1.0))
@test_nowarn sol = solve(prob,Tsit5())
