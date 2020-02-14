using Test

using DiffEqBase
using DiffEqBase: prob2dtmin

@testset "tspan2dtmin" begin
  tspan2dtmin(tspan) = prob2dtmin(ODEProblem(identity, 1, tspan))
  @test tspan2dtmin((10, 100.0)) === eps(100.0)
  @test tspan2dtmin((-10000.0, 100.0)) === eps(10000.0)
  @test tspan2dtmin((1, 2)) === 0
  @test tspan2dtmin((1//10, 2//10)) === 1//10^10
  @test tspan2dtmin((2//10, Inf)) === eps(0.2)
  @test_throws ArgumentError tspan2dtmin((Inf, 100.0))
end
