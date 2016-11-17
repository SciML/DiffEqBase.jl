using DiffEqBase
using Base.Test

@time @testset "Number of Parameters Calculation" begin include("numparameters_test.jl") end
@time @testset "Solution Interface" begin include("solution_get_tests.jl") end
