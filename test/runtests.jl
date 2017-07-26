using DiffEqBase
using Base.Test

@time @testset "Number of Parameters Calculation" begin include("numargs_test.jl") end
@time @testset "Solution Interface" begin include("solution_get_tests.jl") end
@time @testset "Extended Functions" begin include("extended_function_tests.jl") end
@time @testset "Callbacks" begin include("callbacks.jl") end
@time @testset "Constructed Parameterized Functions" begin include("constructed_pf_test.jl") end
@time @testset "Plot Variables" begin include("plot_vars.jl") end
@time @testset "Muladd Macro" begin include("muladd.jl") end
