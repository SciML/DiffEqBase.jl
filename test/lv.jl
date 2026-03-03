using OrdinaryDiffEqTsit5, OrdinaryDiffEqCore, SciMLBase, LinearAlgebra

# Lotka-Volterra equations
function lotka_volterra!(du, u, p, t)
    α, β, δ, γ = 1.5, 1.0, 3.0, 1.0
    x, y = u
    du[1] = α * x - β * x * y
    du[2] = δ * x * y - γ * y
    return nothing
end

# coeffs1 = [0.6825223495861318, -0.4295322984152052]
# coeffs2 = [1.7358772252665537, -1.0070061675696311]

coeffs1 = [2.922772251297381, -2.8028553839288595]
u0 = [1.0, 1.0]
tspan = (0.0, 20.0)
tspan = (0.0, 0.03)

# Record initial signs
initial_conditions = [dot(coeffs1, u0)]
initial_signs = sign.(initial_conditions)

# VCC condition: two linear functions of state
function vcc_condition!(out, u, t, integrator)
    out[1] = dot(coeffs1, u)
    return nothing
end

function vcc_affect!(integrator, event_index)
    @show event_index, integrator.t
    u = integrator.u
    @show integrator.t
    if event_index == 1
        @show integrator.u
        println("Condition value at crossing: ", [dot(coeffs1, u)])
        terminate!(integrator)
    end
    return nothing
end

cb = VectorContinuousCallback(
    vcc_condition!,
    vcc_affect!,
    1;
    rootfind=SciMLBase.RightRootFind,
)

println("Initial conditions: ", initial_conditions)

prob = ODEProblem(lotka_volterra!, u0, tspan)
sol = solve(prob, Tsit5(); callback=cb, abstol=1e-10, reltol=1e-10, dense=true)
sol_u = solve(prob, Tsit5(); abstol=1e-10, reltol=1e-10, dense=true)
# sol = solve(prob, Tsit5(); abstol=1e-10, reltol=1e-10)
sol
nothing

shift(τ, i) =
    if iszero(i)
        τ
    elseif i > 0
        shift(nextfloat(τ), i - 1)
    else
        shift(prevfloat(τ), i + 1)
    end

# 0.23620973794890948

cond2(u) = dot(coeffs2, u)