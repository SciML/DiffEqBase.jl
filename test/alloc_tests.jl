using AllocCheck
using DiffEqBase
using Test

@testset "AllocCheck - Zero Allocations" begin
    @testset "ODE_DEFAULT_NORM" begin
        # Vector{Float64}
        @check_allocs function test_norm_vec_f64(u::Vector{Float64}, t::Float64)
            return DiffEqBase.ODE_DEFAULT_NORM(u, t)
        end
        u = rand(100)
        t = 0.0
        @test test_norm_vec_f64(u, t) isa Float64

        # Scalar Float64
        @check_allocs function test_norm_scalar_f64(u::Float64, t::Float64)
            return DiffEqBase.ODE_DEFAULT_NORM(u, t)
        end
        @test test_norm_scalar_f64(1.5, 0.0) isa Float64

        # Vector{ComplexF64}
        @check_allocs function test_norm_vec_c64(u::Vector{ComplexF64}, t::Float64)
            return DiffEqBase.ODE_DEFAULT_NORM(u, t)
        end
        u_complex = rand(ComplexF64, 100)
        @test test_norm_vec_c64(u_complex, 0.0) isa Float64

        # Scalar ComplexF64
        @check_allocs function test_norm_scalar_c64(u::ComplexF64, t::Float64)
            return DiffEqBase.ODE_DEFAULT_NORM(u, t)
        end
        @test test_norm_scalar_c64(1.5 + 0.5im, 0.0) isa Float64
    end

    @testset "calculate_residuals!" begin
        # In-place version with Vector{Float64}
        @check_allocs function test_residuals_inplace_f64(
            out::Vector{Float64},
            ũ::Vector{Float64},
            u₀::Vector{Float64},
            u₁::Vector{Float64},
            α::Float64,
            ρ::Float64,
            t::Float64
        )
            return DiffEqBase.calculate_residuals!(
                out, ũ, u₀, u₁, α, ρ,
                DiffEqBase.ODE_DEFAULT_NORM, t, DiffEqBase.False()
            )
        end

        n = 100
        out = zeros(n)
        ũ = rand(n)
        u₀ = rand(n)
        u₁ = rand(n)
        @test test_residuals_inplace_f64(out, ũ, u₀, u₁, 1e-6, 1e-3, 0.0) === nothing
    end

    @testset "calculate_residuals (scalar)" begin
        @check_allocs function test_residuals_scalar_f64(
            ũ::Float64,
            u₀::Float64,
            u₁::Float64,
            α::Float64,
            ρ::Float64,
            t::Float64
        )
            return DiffEqBase.calculate_residuals(
                ũ, u₀, u₁, α, ρ,
                DiffEqBase.ODE_DEFAULT_NORM, t
            )
        end
        @test test_residuals_scalar_f64(0.001, 1.0, 1.1, 1e-6, 1e-3, 0.0) isa Float64
    end

    @testset "UNITLESS_ABS2" begin
        @check_allocs function test_unitless_abs2_f64(u::Vector{Float64})
            return DiffEqBase.UNITLESS_ABS2(u)
        end
        u = rand(100)
        @test test_unitless_abs2_f64(u) isa Float64

        # Scalar
        @check_allocs function test_unitless_abs2_scalar(u::Float64)
            return DiffEqBase.UNITLESS_ABS2(u)
        end
        @test test_unitless_abs2_scalar(1.5) isa Float64
    end

    @testset "recursive_length" begin
        @check_allocs function test_recursive_length_f64(u::Vector{Float64})
            return DiffEqBase.recursive_length(u)
        end
        u = rand(100)
        @test test_recursive_length_f64(u) == 100
    end

    @testset "NAN_CHECK" begin
        @check_allocs function test_nan_check_f64(u::Vector{Float64})
            return DiffEqBase.NAN_CHECK(u)
        end
        u = rand(100)
        @test test_nan_check_f64(u) == false

        u_nan = rand(100)
        u_nan[50] = NaN
        @test test_nan_check_f64(u_nan) == true

        # Scalar
        @check_allocs function test_nan_check_scalar(u::Float64)
            return DiffEqBase.NAN_CHECK(u)
        end
        @test test_nan_check_scalar(1.5) == false
        @test test_nan_check_scalar(NaN) == true
    end

    @testset "INFINITE_OR_GIANT" begin
        @check_allocs function test_infinite_check_f64(u::Vector{Float64})
            return DiffEqBase.INFINITE_OR_GIANT(u)
        end
        u = rand(100)
        @test test_infinite_check_f64(u) == false

        u_inf = rand(100)
        u_inf[50] = Inf
        @test test_infinite_check_f64(u_inf) == true

        # Scalar
        @check_allocs function test_infinite_check_scalar(u::Float64)
            return DiffEqBase.INFINITE_OR_GIANT(u)
        end
        @test test_infinite_check_scalar(1.5) == false
        @test test_infinite_check_scalar(Inf) == true
    end

    @testset "findall_events!" begin
        # In-place with Vector{Float64}
        affect! = (integrator) -> nothing
        affect_neg! = (integrator) -> nothing

        @check_allocs function test_findall_events_f64(
            next_sign::Vector{Float64},
            prev_sign::Vector{Float64}
        )
            return DiffEqBase.findall_events!(
                next_sign,
                affect!,
                affect_neg!,
                prev_sign
            )
        end

        n = 100
        next_sign = randn(n)
        prev_sign = randn(n)
        @test test_findall_events_f64(next_sign, prev_sign) === next_sign
    end
end
