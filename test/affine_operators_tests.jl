using DiffEqBase
using Test
using Random

mutable struct TestDiffEqOperator{T} <: DiffEqBase.AbstractDiffEqLinearOperator{T}
    m::Int
    n::Int
end

TestDiffEqOperator(A::AbstractMatrix{T}) where {T} =
    TestDiffEqOperator{T}(size(A)...)

Base.size(A::TestDiffEqOperator) = (A.m, A.n)


A = TestDiffEqOperator([0 0; 0 1])
B = TestDiffEqOperator([0 0 0; 0 1 0; 0 0 2])

@test_throws ErrorException AffineDiffEqOperator{Int64}((A,B),())

@testset "DiffEq linear operators" begin
    Random.seed!(0)
    M = rand(2,2); A = DiffEqArrayOperator(M)
    b = rand(2)
    u = rand(2)
    p = rand(1)
    t = rand()
    As_list = [(A,), (A, A)]#, (A, α)]
    bs_list = [(), (b,), (2b,), (b, 2b)]
    @testset "combinations of A's and b's" for As in As_list, bs in bs_list
        L = AffineDiffEqOperator{Float64}(As, bs, zeros(2))
        mysum = sum(A*u for A in As)
        for b in bs; mysum .+= b; end
        @test L(u,p,t) == mysum
    end
end
