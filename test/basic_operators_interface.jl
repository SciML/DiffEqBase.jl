using DiffEqBase, Random, LinearAlgebra, Test

@testset "Identity Operators" begin
  u = [1.0, 2.0]; du = zeros(2)
  Id = DiffEqIdentity(u)
  @test Id * u == u
  mul!(du, Id, u); @test du == u
  @test size(Id) == (2,2)
end

@testset "Scalar Operators" begin
  u = [1.0, 2.0]; u2 = [1.0, 2.0]
  α = DiffEqScalar(2.0)
  @test convert(Number, α) == 2.0
  @test α * u == 2.0u
  lmul!(α, u2); @test u2 == 2.0u
  @test size(α) == ()
  @test is_constant(α) == true
end

@testset "Array Operators" begin
  Random.seed!(0); A = rand(2,2); u = rand(2); du = zeros(2)
  L = DiffEqArrayOperator(A)
  @test Matrix(L) == A
  @test size(L) == size(A)
  @test L * u == A * u
  mul!(du, L, u); @test du == A*u
  @test lu(L) \ u ≈ A \ u
  @test opnorm(L) == opnorm(A)
  @test exp(L) == exp(A)
  @test L[1,2] == A[1,2]
  @test is_constant(L) == true
end

@testset "Mutable Array Operators" begin
  Random.seed!(0); A = rand(2,2); u = rand(2); du = zeros(2)
  update_func = (_A,u,p,t) -> _A .= t * A
  Lt = DiffEqArrayOperator(zeros(2,2); update_func=update_func)
  t = 5.0
  @test is_constant(Lt) == false
  @test Lt(u,nothing,t) ≈ (t*A) * u
  Lt(du,u,nothing,t); @test du ≈ (t*A) * u
end
