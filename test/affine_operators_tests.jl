using DiffEqBase
using Test
using Random

mutable struct TestSciMLOperator{T} <: DiffEqBase.AbstractSciMLLinearOperator{T}
    m::Int
    n::Int
end

TestSciMLOperator(A::AbstractMatrix{T}) where {T} =
    TestSciMLOperator{T}(size(A)...)

Base.size(A::TestSciMLOperator) = (A.m, A.n)


A = TestSciMLOperator([0 0; 0 1])
B = TestSciMLOperator([0 0 0; 0 1 0; 0 0 2])

@test_throws ErrorException AffineSciMLOperator{Int64}((A,B),())

@testset "DiffEq linear operators" begin
    Random.seed!(0)
    M = rand(2,2); A = MatrixOperator(M)
    b = rand(2)
    u = rand(2)
    p = rand(1)
    t = rand()
    As_list = [(A,), (A, A)]#, (A, α)]
    bs_list = [(), (b,), (2b,), (b, 2b)]
    @testset "combinations of A's and b's" for As in As_list, bs in bs_list
        L = AffineSciMLOperator{Float64}(As, bs, zeros(2))
        mysum = sum(A*u for A in As)
        for b in bs; mysum .+= b; end
        @test L(u,p,t) == mysum
    end
end
