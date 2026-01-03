# Precompilation workload for DiffEqBase
# This precompiles commonly used code paths to reduce TTFX

PrecompileTools.@setup_workload begin
    # Minimal setup for precompilation
    u0 = [1.0, 0.0, 0.0]
    u1 = [1.1, 0.1, 0.1]
    p = [1.0, 2.0, 3.0]
    t = 0.0
    α = 1e-6
    ρ = 1e-3

    PrecompileTools.@compile_workload begin
        # Precompile ODE_DEFAULT_NORM for Vector{Float64} (most common case)
        ODE_DEFAULT_NORM(u0, t)

        # Precompile NAN_CHECK for Vector{Float64}
        NAN_CHECK(u0)

        # Precompile INFINITE_OR_GIANT for Vector{Float64}
        INFINITE_OR_GIANT(u0)

        # Precompile calculate_residuals for Vector{Float64} (most common case)
        calculate_residuals(u0, u1, α, ρ, ODE_DEFAULT_NORM, t)

        # Precompile in-place version
        out = similar(u0)
        calculate_residuals!(out, u0, u1, α, ρ, ODE_DEFAULT_NORM, t)

        # Precompile with explicit differences
        ũ = u1 .- u0
        calculate_residuals(ũ, u0, u1, α, ρ, ODE_DEFAULT_NORM, t)
        calculate_residuals!(out, ũ, u0, u1, α, ρ, ODE_DEFAULT_NORM, t)
    end
end
