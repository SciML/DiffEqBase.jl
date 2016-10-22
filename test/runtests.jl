using DiffEqBase
using Base.Test

# write your own tests here
println("Number of Parameters Calc Tests")
@time @test include("numparameters_test.jl")
println("Solver Interface Tests")
@time @test include("solution_get_tests.jl")
