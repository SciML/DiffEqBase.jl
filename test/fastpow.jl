using DiffEqBase: fastlog2, fastpow
using Test

@testset "Fast log2" begin
    for x in 0.001:0.001:1.2 # (0, 1+something] is the domain which a controller uses
        @test log2(x)≈fastlog2(Float32(x)) atol=1e-3
    end
end

@testset "Fast pow" begin
    @test fastpow(1, 1) isa Float64
    @test fastpow(1.0, 1.0) isa Float64
    errors = [abs(^(x, y) - fastpow(x, y)) for x in 0.001:0.001:1, y in 0.08:0.001:0.5]
    @test maximum(errors) < 1e-4
end