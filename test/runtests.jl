using SafeTestsets

const GROUP = get(ENV, "GROUP", "All")
const is_APPVEYOR = ( Sys.iswindows() && haskey(ENV,"APPVEYOR") )
const is_TRAVIS = haskey(ENV,"TRAVIS")

@time begin
if GROUP == "All" || GROUP == "Core"
    @time @safetestset "Fast Power" begin include("fastpow.jl") end
    @time @safetestset "Fast Broadcast" begin include("fastbc.jl") end
    @time @safetestset "Number of Parameters Calculation" begin include("numargs_test.jl") end
    @time @safetestset "Data Arrays" begin include("data_array_tests.jl") end
    @time @safetestset "Existence functions" begin include("existence_functions.jl") end
    @time @safetestset "Callbacks" begin include("callbacks.jl") end
    @time @safetestset "Plot Variables" begin include("plot_vars.jl") end
    @time @safetestset "Problem Creation Tests" begin include("problem_creation_tests.jl") end
    @time @safetestset "Remake tests" begin include("remake_tests.jl") end
    @time @safetestset "Affine differential equation operators" begin include("affine_operators_tests.jl") end
    @time @safetestset "TableTraits" begin include("tabletraits_tests.jl") end
    @time @safetestset "Integrator interface" begin include("integrator_tests.jl") end
    @time @safetestset "Export tests" begin include("export_tests.jl") end
    @time @safetestset "High Level solve Interface" begin include("high_level_solve.jl") end
    @time @safetestset "DiffEqFunction tests" begin include("diffeqfunction_tests.jl") end
    @time @safetestset "Internal Euler" begin include("internal_euler_test.jl") end
    @time @safetestset "Numargs" begin include("numargs_test.jl") end
    @time @safetestset "Basic Operators Interface" begin include("basic_operators_interface.jl") end
    @time @safetestset "Norm" begin include("norm.jl") end
    @time @safetestset "Utils" begin include("utils.jl") end
end

if !is_APPVEYOR && GROUP == "Downstream"
    # add additional packages
    using Pkg
    Pkg.activate("downstream")
    Pkg.develop(PackageSpec(path=joinpath(pwd(), "..")))
    Pkg.instantiate()

    @time @safetestset "Null Parameters" begin include("downstream/null_params_test.jl") end
    @time @safetestset "Ensemble Simulations" begin include("downstream/ensemble.jl") end
    @time @safetestset "Ensemble Analysis" begin include("downstream/ensemble_analysis.jl") end
    @time @safetestset "Default linsolve with structure" begin include("downstream/default_linsolve_structure.jl") end
    @time @safetestset "Callback Merging Tests" begin include("downstream/callback_merging.jl") end
    @time @safetestset "LabelledArrays Tests" begin include("downstream/labelledarrays.jl") end
    @time @safetestset "ODE Event Tests" begin include("downstream/ode_event_tests.jl") end
    @time @safetestset "Event Detection Tests" begin include("downstream/event_detection_tests.jl") end
    @time @safetestset "PSOS and Energy Conservation Event Detection" begin include("downstream/psos_and_energy_conservation.jl") end
    @time @safetestset "DE stats" begin include("downstream/destats_tests.jl") end
    @time @safetestset "DEDataArray" begin include("downstream/data_array_regression_tests.jl") end
    @time @safetestset "Concrete_solve Tests" begin include("downstream/concrete_solve_tests.jl") end
end

if !is_APPVEYOR && GROUP == "GPU"
  @time @safetestset "Simple GPU" begin include("gpu/simple_gpu.jl") end
end

end
