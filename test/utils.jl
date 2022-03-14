using Test

using DiffEqBase
using DiffEqBase: prob2dtmin, timedepentdtmin

@testset "tspan2dtmin" begin
  tspan2dtmin(tspan; kwargs...) = prob2dtmin(ODEProblem(identity, 1, tspan); kwargs...)
  @test tspan2dtmin((10, 100.0)) === eps(100.0)
  @test tspan2dtmin((-10000.0, 100.0)) === eps(10000.0)
  @test tspan2dtmin((1, 2)) === 0
  @test tspan2dtmin((1//10, 2//10)) === 1//10^10
  @test tspan2dtmin((2//10, Inf)) === eps(1.0)
  @test tspan2dtmin((2//1, Inf)) === eps(2.0)
  @test tspan2dtmin((0, Inf)) === tspan2dtmin((0.0, Inf)) === eps(1.0)
  @test_throws ArgumentError tspan2dtmin((Inf, 100.0))
  @test tspan2dtmin((0f0, 1f5); use_end_time=false) === eps(1f0)
  @test timedepentdtmin(10f0, eps(1f0)) === eps(10f0)
  @test timedepentdtmin(10, eps(1f0)) === eps(1f0)
end

@testset "prob2dtmin" begin
  @test prob2dtmin((0.0, 10.0), 1.0, false) == eps(Float64)
  @test prob2dtmin((0f0, 10f0), 1f0, false) == eps(Float32)
  @test prob2dtmin((0.0, 10.0), ForwardDiff.Dual(1.0), false) == eps(Float64)
end
