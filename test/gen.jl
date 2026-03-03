using OrdinaryDiffEqTsit5, OrdinaryDiffEqCore, SciMLBase, LinearAlgebra

# Lotka-Volterra equations
function lotka_volterra!(du, u, p, t)
    α, β, δ, γ = 1.5, 1.0, 3.0, 1.0
    x, y = u
    du[1] = α * x - β * x * y
    du[2] = δ * x * y - γ * y
    return nothing
end

function run_once(; seed=nothing)
    if seed !== nothing
        Random.seed!(seed)
    end

    # Random coefficients for two linear conditions: c' * u
    coeffs1 = randn(2)

    u0 = [1.0, 1.0]
    tspan = (0.0, 20.0)

    # Record initial signs
    initial_signs = [sign(dot(coeffs1, u0))]

    # VCC condition: two linear functions of state
    function vcc_condition!(out, u, t, integrator)
        out[1] = dot(coeffs1, u)
        return nothing
    end

    function vcc_affect!(integrator, event_index)
        u = integrator.u
        vals = [dot(coeffs1, u)]
        v = vals[event_index]
        if !iszero(v) && sign(v) == initial_signs[event_index]
            @show coeffs1 u event_index v initial_signs
            error("VCC fired but value has same sign as initial — RightRootFind bug?")
        else
            # termine simulation
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

    prob = ODEProblem(lotka_volterra!, u0, tspan)
    sol = solve(prob, Tsit5(); callback=cb, abstol=1e-10, reltol=1e-10)
    return sol
end

# Main loop — fish for the bug
i = 0
while true
    global i += 1
    if i % 1000 == 0
        println("Iteration $i ...")
    end
    try
        run_once()
    catch e
        if e isa ErrorException && contains(e.msg, "RightRootFind")
            println("\n*** Bug found at iteration $i ***")
            rethrow()
        else
            rethrow()
        end
    end
end

println("No bug found after 100_000 iterations.")
