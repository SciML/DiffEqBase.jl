using DiffEqBase
using Base.Test

tic()
@time @testset "Number of Parameters Calculation" begin include("numargs_test.jl") end
@time @testset "Data Arrays" begin include("data_array_tests.jl") end
@time @testset "Extended Functions" begin include("extended_function_tests.jl") end
@time @testset "DiffEqFunctions" begin include("diffeqfunction_tests.jl") end
@time @testset "Callbacks" begin include("callbacks.jl") end
@time @testset "Plot Variables" begin include("plot_vars.jl") end
@time @testset "Problem Creation Tests" begin include("problem_creation_tests.jl") end
@time @testset "Affine differential equation operators" begin include("affine_operators_tests.jl") end
toc()
