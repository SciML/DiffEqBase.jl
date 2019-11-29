# https://github.com/JuliaDiffEq/DifferentialEquations.jl/issues/525

using OrdinaryDiffEq, StaticArrays

mutable struct SimType{T} <: DEDataVector{T}
    x::Array{T,1}
    f1::T
end

function f(u,p,t) # new out-of-place definition
    SimType([-0.5*u[1] + u.f1,
    -0.5*u[2]],u.f1)
end

function f!(du,u,p,t) # old in-place definition
    du[1] = -0.5*u[1] + u.f1
    du[2] = -0.5*u[2]
end

const tstop1 = [5.]
const tstop2 = [8.]

function condition(u,t,integrator)
  t in tstop1
end

function condition2(u,t,integrator)
  t in tstop2
end

function affect!(integrator)
  for c in full_cache(integrator)
    c.f1 = 1.5
  end
end

function affect2!(integrator)
  for c in full_cache(integrator)
    c.f1 = -1.5
  end
end

save_positions = (true,true)
cb = DiscreteCallback(condition, affect!, save_positions=save_positions)
save_positions = (false,true)
cb2 = DiscreteCallback(condition2, affect2!, save_positions=save_positions)
cbs = CallbackSet(cb,cb2)
u0 = SimType([10.0;10.0], 0.0)

prob_inplace = ODEProblem(f!,u0,(0.0,10.0))
prob = ODEProblem(f,u0,(0.0,10.0))

const tstop = [5.;8.]

sol = solve(prob_inplace,Tsit5(),callback = cbs, tstops=tstop)
sol = solve(prob,Tsit5(),callback = cbs, tstops=tstop)

# https://github.com/JuliaDiffEq/DifferentialEquations.jl/issues/336

const A = SMatrix{2,2}([0 1;
                        0 0])

mutable struct MyStruct{T} <: DEDataVector{T}
    x::MVector{2,T}
    a::SVector{2,T}
end

function dyn(du,u,t,p)
    u.a = SVector{2}(0.0, 0.1)
    du .= A*u.x + u.a

    return nothing
end

u0   = MyStruct(MVector{2}(0.,0.), SVector{2}(0.,0.))
prob = ODEProblem(dyn, u0, (0.,10.))

@test copy(u0) isa MyStruct
@test zero(u0) isa MyStruct
@test similar(u0) isa MyStruct
@test similar(u0,Float64) isa MyStruct
@test similar(u0,Float64,size(u0)) isa MyStruct

sol  = solve(prob, Tsit5())


# https://github.com/JuliaDiffEq/StochasticDiffEq.jl/issues/247

using OrdinaryDiffEq
using StochasticDiffEq

mutable struct SimType{T} <: DEDataVector{T}
    x::Array{T,1}
    f1::T
end

function f(du,u,p,t)
    du[1] = -0.5*u[1] + u.f1
    du[2] = -0.5*u[2]
end

const tstop1 = [5.]
const tstop2 = [8.]


function condition(u,t,integrator)
  t in tstop1
end

function condition2(u,t,integrator)
  t in tstop2
end

function affect!(integrator)
  for c in full_cache(integrator)
    c.f1 = 1.5
  end
end

function affect2!(integrator)
  for c in full_cache(integrator)
    c.f1 = -1.5
  end
end

save_positions = (true,true)

cb = DiscreteCallback(condition, affect!, save_positions=save_positions)

save_positions = (false,true)

cb2 = DiscreteCallback(condition2, affect2!, save_positions=save_positions)

cbs = CallbackSet(cb,cb2)

u0 = SimType([10.0;10.0], 0.0)
prob = ODEProblem(f,u0,(0.0,10.0))

const tstop = [5.;8.]

# here the new part

function g(du,u,p,t)
  du[1] = 1.0
  du[2] = 1.2
end

dt = 1/2^4

prob2 = SDEProblem(f,g,u0,(0.0,10.0))

# this creates an error

sol = solve(prob2,callback = cbs,tstops=tstop,EM(),dt=dt,saveat=collect(8:0.1:10))

# https://github.com/JuliaDiffEq/DiffEqBase.jl/issues/327

struct SimulationState{T,U} <: DiffEqBase.DEDataArray{T, 1}
    x::Vector{T}
    u::U
end

#function Base.convert(::Type{SimulationState{T, U}}, a::AbstractArray{T}) where {T,U}
#    SimulationState(a, zero(eltype(a)))
#end

function open_loop(state, parameters, time)
    v = state[1] * -0.1 + state.u
    return SimulationState([v],state.u)
end

initial_conditions = SimulationState([10.0], 0.0)
time_span = (0.0, 20.0)
ode_prob = ODEProblem(open_loop, initial_conditions, time_span, nothing)

sol = solve(ode_prob, Tsit5(), reltol=1e-8, abstol=1e-8)
