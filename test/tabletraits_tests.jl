using DiffEqBase, DiffEqBase.InternalEuler
using IteratorInterfaceExtensions
using Base.Test

f_1dlinear = (u,p,t) -> 1.01u
prob = ODEProblem(f_1dlinear,rand(),(0.0,1.0))
sol =solve(prob,InternalEuler.FwdEulerAlg();dt=1//2^(4))
array = collect(getiterator(sol))

@test length(array) == 17
@test collect(fieldnames(array[1])) == [:timestamp, :value]

f_2dlinear = (du,u,p,t) -> du.=1.01u
prob = ODEProblem(f_2dlinear,rand(2),(0.0,1.0))
sol =solve(prob,InternalEuler.FwdEulerAlg();dt=1//2^(4))
array = collect(getiterator(sol))

@test length(array) == 17
@test collect(fieldnames(array[1])) == [:timestamp, :value1, :value2]

struct LotkaVolterra
    syms::Vector{Symbol}
end
function (::LotkaVolterra)(du,u,p,t)
    du[1] = 1.5*u[1] - 1*u[1]*u[2]
    du[2] = -3*u[2] + 1*u[1]*u[2]
end
f_2dlinear_named = LotkaVolterra([:x,:y])
prob = ODEProblem(f_2dlinear_named,[1.0,1.0],(0.0,1.0))
sol =solve(prob,InternalEuler.FwdEulerAlg(),dt=0.1)
array = collect(getiterator(sol))

@test length(array) == 11
@test collect(fieldnames(array[1])) == [:timestamp, :x, :y]
