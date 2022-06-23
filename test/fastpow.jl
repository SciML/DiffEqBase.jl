using DiffEqBase: fastlog2, _exp2, fastpow
using Test

@testset "Fast log2" begin
    for x = 0.001:0.001:1.2 # (0, 1+something] is the domain which a controller uses
        @test log2(x) ≈ fastlog2(Float32(x)) atol = 1e-3
    end
end

@testset "Exp2" begin
    for x = -100:0.01:3
        @test exp2(x) ≈ _exp2(Float32(x)) atol = 1e-6
    end
end

@testset "Fast pow" begin
    @test fastpow(1, 1) isa Float32
    @test fastpow(1.0, 1.0) isa Float32
    errors = [abs(^(x, y) - fastpow(x, y)) for x = 0.001:0.001:1, y = 0.08:0.001:0.5]
    @test maximum(errors) < 1e-4
end
