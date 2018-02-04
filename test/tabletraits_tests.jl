using DiffEqBase, DiffEqBase.InternalEuler
using ParameterizedFunctions
using IteratorInterfaceExtensions
using Base.Test

f_1dlinear = (u,p,t) -> 1.01u
prob = ODEProblem(f_1dlinear,rand(),(0.0,1.0))
sol =solve(prob,InternalEuler.FwdEulerAlg();dt=1//2^(4))
array = collect(getiterator(sol))

@test length(array) == 17
@test fieldnames(array[1]) == [:timestamp, :value]

f_2dlinear = (du,u,p,t) -> du.=1.01u
prob = ODEProblem(f_2dlinear,rand(2),(0.0,1.0))
sol =solve(prob,InternalEuler.FwdEulerAlg();dt=1//2^(4))
array = collect(getiterator(sol))

@test length(array) == 17
@test fieldnames(array[1]) == [:timestamp, :value1, :value2]

f_2dlinear_named = @ode_def LotkaVolterra begin
    dx = 1.5*x - 1*x*y
    dy = -3*y + 1*x*y
end
  
prob = ODEProblem(f_2dlinear_named,[1.0,1.0],(0.0,1.0))
sol =solve(prob,Tsit5())
array = collect(getiterator(sol))

@test length(array) == 7
@test fieldnames(array[1]) == [:timestamp, :x, :y]
