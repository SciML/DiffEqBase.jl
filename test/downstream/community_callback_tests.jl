using OrdinaryDiffEq, DiffEqCallbacks, LinearAlgebra

# https://github.com/SciML/DiffEqBase.jl/issues/564 : Fixed
gravity = 9.8
stiffness = 500
equilibrium_length = 1
T = 5.0

f(u, p, t) = begin
    x1, x2, dx1, dx2 = u
    length = abs(x2 - x1)
    spring_force = stiffness * (equilibrium_length - length)
    ddx1 = -gravity - spring_force
    ddx2 = -gravity + spring_force
    if x1 <= 0
        ddx1 = max(0, ddx1)
    end
    if x2 <= 0
        ddx2 = max(0, ddx2)
    end
    [dx1, dx2, ddx1, ddx2]
end

sol = solve(
    ODEProblem(f, [5.0, 6.0, 0.0, 0.0], (0.0, T)),
    # Euler(),
    # dt=0.005,
    Rosenbrock23(),
    callback = ContinuousCallback(
        (u, _, _) -> u[1],
        (integrator) -> (integrator.u[1] = 0; integrator.u[3] = 0),
    ),
    # callback = ContinuousCallback((u, _, _) -> u[1], (integrator) -> (integrator.u[3] = 0)),
    reltol = 1e-3,
    abstol = 1e-3,
)

@show sol.destats

# https://github.com/SciML/DiffEqBase.jl/issues/553 : Floating point issue is resolved but some other error occurss
function model(du, u, p, t)
    du[1] = 0.0
    for i = 2:(length(du)-1)
        du[i] = p[i] * (u[i-1] - u[i])
    end
    du[end] = p[end] * (p[1] * u[end-1] - u[end])
    return nothing
end

perror = [
    1.0,
    0.02222434508140991,
    0.017030281542289794,
    0.015917011145559996,
    0.1608874463597176,
    0.13128016561792297,
    0.11056834258380167,
    0.5222141958458832,
    1.0711942201995688,
    0.2672878398678257,
    8.900058706990183,
    0.010760065201065117,
    0.016319181296867765,
    2.2693845639611925,
    0.2152216345154439,
    0.029186712540925457,
    0.21419429135100806,
    0.029177617589788596,
    0.03064986043089549,
    0.023280222517122397,
    6.931251277770224,
]
y_max = 0.002604806609572015
u0 = [1, zeros(length(perror) - 1)...]
tspan = (0.0, 5000.0)

condition(u, t, i) = (t == 1.0)
affect!(i) = (i.u[1] = 0.0)

condition2(u, t, i) = u[end] - y_max / 2.0
t_half_1 = 0.0
affect2!(i) = (t_half_1 = i.t)

prob = ODEProblem(model, u0, tspan, perror)
sol = solve(
    prob,
    Rosenbrock23();
    callback = CallbackSet(
        PositiveDomain(),
        DiscreteCallback(condition, affect!),
        ContinuousCallback(condition2, affect2!, terminate!),
    ),
    tstops = [1.0],
    force_dtmin = true,
)

# https://github.com/SciML/DiffEqBase.jl/issues/515 : Fixed

using StaticArrays
using MultiScaleArrays

t_last = 0.0
function attactor(du, u, p, t)
    α, β = p
    n = length(u.nodes)
    return for k = 1:n
        du.nodes[k] = zero(du.nodes[k])
        for j = 1:n
            if (k == j)
                du.nodes[k] .+=
                    [u.nodes[k][3], u.nodes[k][4], -β * u.nodes[k][3], -β * u.nodes[k][4]]
            else
                du.nodes[k][3:4] .+= α * (u.nodes[j][1:2] - u.nodes[k][1:2])
            end
        end
    end
end

struct Thingy{B} <: AbstractMultiScaleArrayLeaf{B}
    values::Vector{B}
end

struct PhysicsLaw{T<:AbstractMultiScaleArray,B<:Number} <: AbstractMultiScaleArrayHead{B}
    nodes::Vector{T}
    values::Vector{B}
    end_idxs::Vector{Int}
end

Newton = construct(
    PhysicsLaw,
    [
        Thingy([-700.0, -350.0, 0.0, 0.0]),
        Thingy([-550.0, -150.0, 0.0, 0.0]),
        Thingy([-600.0, 15.0, 0.0, 10.0]),
        Thingy([200.0, -200.0, 5.0, -5.0]),
    ],
)

parameters = [1e-2, 0.06]


function condition(out, u, t, integrator)
    i = 0
    n = length(u.nodes)
    for k = 1:n
        for l = (k+1):n
            i += 1
            out[i] = sum(abs2, u.nodes[k][1:2] .- u.nodes[l][1:2]) - 10000
        end
    end
end

function affect!(integrator, idx)
    i = 0
    u = integrator.u
    n = length(u.nodes)
    return for k = 1:n
        for l = (k+1):n
            i += 1
            if idx == i
                x₁ = u.nodes[k][1:2]
                v₁ = u.nodes[k][3:4]
                x₂ = u.nodes[l][1:2]
                v₂ = u.nodes[l][3:4]
                # https://stackoverflow.com/a/35212639
                v₁ = (
                    v₁ -
                    2 / (1 + 1) * (dot(v₁ - v₂, x₁ - x₂) / sum(abs2, x₁ - x₂) * (x₁ - x₂))
                )
                v₂ = -(
                    v₂ -
                    2 / (1 + 1) * (dot(v₂ - v₁, x₂ - x₁) / sum(abs2, x₂ - x₁) * (x₂ - x₁))
                )

                println("Collision handeled.")

                m = (x₁ + x₂) / 2


                u.nodes[k][3:4] .= v₁
                u.nodes[l][3:4] .= v₂

                set_u!(integrator, u)
                println(sqrt(sum(abs2, x₁ .- x₂)) - 100, ":", v₁ ./ v₂)
                println(
                    norm(v₁),
                    ":",
                    norm(v₂),
                    ":",
                    integrator.t,
                    ":",
                    integrator.t - t_last,
                )
                global t_last = integrator.t
                break
            end
        end
    end
end

cback = VectorContinuousCallback(
    condition,
    affect!,
    (x -> Int(((x - 1) * x) / 2))(length(Newton.nodes)), #Fix  was here
)


problemp = ODEProblem(attactor, Newton, (0.0, Inf), parameters)

world = init(problemp, AutoTsit5(Rosenbrock23()); save_everystep = false, callback = cback)

dt = 0.2

for i = 1:1000
    step!(world, dt)
end

## https://github.com/SciML/OrdinaryDiffEq.jl/issues/1528

function f!(out, u, p, t)
    out[1] = 0
    out[2] = u[3]
    out[3] = -1.0 * (u[2] - u[1])
end
u0 = [0, 0, 1.0]
function cond!(out, u, t, i)
    out[1] = u[3]
    nothing
end
function affect!(int, idx)
    terminate!(int)
end
cb = VectorContinuousCallback(cond!, affect!, nothing, 1)

u0 = [0.0, 0.0, 1.0]
prob = ODEProblem(f!, u0, (0.0, 10.0); callback = cb)
soln = solve(prob, Tsit5())
@test soln.t[end] ≈ 4.712347213360699
