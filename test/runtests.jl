using SafeTestsets

if haskey(ENV,"GROUP")
    group = ENV["GROUP"]
else
    group = "All"
end

is_APPVEYOR = ( Sys.iswindows() && haskey(ENV,"APPVEYOR") )
is_TRAVIS = haskey(ENV,"TRAVIS")

@time begin
if group == "All" || group == "Core"
    @time @safetestset "Fast Broadcast" begin include("fastbc.jl") end
    @time @safetestset "Number of Parameters Calculation" begin include("numargs_test.jl") end
    @time @safetestset "Data Arrays" begin include("data_array_tests.jl") end
    @time @safetestset "Extended Functions" begin include("extended_function_tests.jl") end
    @time @safetestset "Callbacks" begin include("callbacks.jl") end
    @time @safetestset "Plot Variables" begin include("plot_vars.jl") end
    @time @safetestset "Problem Creation Tests" begin include("problem_creation_tests.jl") end
    @time @safetestset "Remake tests" begin include("remake_tests.jl") end
    @time @safetestset "Affine differential equation operators" begin include("affine_operators_tests.jl") end
    @time @safetestset "TableTraits" begin include("tabletraits_tests.jl") end
    @time @safetestset "Integrator interface" begin include("integrator_tests.jl") end
    @time @safetestset "Export tests" begin include("export_tests.jl") end
    @time @safetestset "High Level solve Interface" begin include("high_level_solve.jl") end
end
if !is_APPVEYOR && group == "Downstream"
    if is_TRAVIS
      using Pkg
      Pkg.add("OrdinaryDiffEq")
      Pkg.add("StochasticDiffEq")
      Pkg.add("DiffEqProblemLibrary")
    end
    @time @safetestset "Ensemble Simulations" begin include("downstream/ensemble.jl") end
    @time @safetestset "Ensemble Analysis" begin include("downstream/ensemble_analysis.jl") end
    @time @safetestset "ODE Event Tests" begin include("downstream/ode_event_tests.jl") end
    @time @safetestset "Event Detection Tests" begin include("downstream/event_detection_tests.jl") end
    @time @safetestset "PSOS and Energy Conservation Event Detection" begin include("downstream/psos_and_energy_conservation.jl") end
end
end
