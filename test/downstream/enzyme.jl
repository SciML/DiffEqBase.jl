using Enzyme, EnzymeTestUtils
using DiffEqBase: fastlog2, fastpow
using Test

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