using DiffEqBase: fastlog2, _exp2, fastpow
using Enzyme, EnzymeTestUtils
using Test

@testset "Fast log2" begin
    for x in 0.001:0.001:1.2 # (0, 1+something] is the domain which a controller uses
        @test log2(x)≈fastlog2(Float32(x)) atol=1e-3
    end
end

@testset "Exp2" begin
    for x in -100:0.01:3
        @test exp2(x)≈_exp2(Float32(x)) atol=1e-6
    end
end

@testset "Fast pow" begin
    @test fastpow(1, 1) isa Float32
    @test fastpow(1.0, 1.0) isa Float32
    errors = [abs(^(x, y) - fastpow(x, y)) for x in 0.001:0.001:1, y in 0.08:0.001:0.5]
    @test maximum(errors) < 1e-4
end

@testset "Fast pow - Enzyme forward rule" begin
    @testset for RT in (Duplicated, DuplicatedNoNeed),
        Tx in (Const, Duplicated),
        Ty in (Const, Duplicated)
        x = 3.0
        y = 2.0
        test_forward(fastpow, RT, (x, Tx), (y, Ty), atol=0.005, rtol=0.005)
    end
end

@testset "Fast pow - Enzyme reverse rule" begin
    @testset for RT in (Active,),
        Tx in (Active,),
        Ty in (Active,)
        x = 2.0
        y = 3.0
        test_reverse(fastpow, RT, (x, Tx), (y, Ty), atol=0.001, rtol=0.001)
    end
end