using DifferentialEquationsBase
using Base.Test

# write your own tests here
println("Number of Parameters Calc Tests")
@time @test include("internals/numparameters_test.jl")
println("Solver Interface Tests")
@time @test include("internals/solution_get_tests.jl")
