using BenchmarkTools, CUDA, DiffEqBase, Test
CUDA.allowscalar(false)

du = cu(rand(4))
u = cu(rand(4))
uprev = cu(rand(4))

const TERMINATION_CONDITIONS = [
    SteadyStateDiffEqTerminationMode(), SimpleNonlinearSolveTerminationMode(),
    NormTerminationMode(), RelTerminationMode(), RelNormTerminationMode(),
    AbsTerminationMode(), AbsNormTerminationMode(), RelSafeTerminationMode(),
    AbsSafeTerminationMode(), RelSafeBestTerminationMode(), AbsSafeBestTerminationMode()
]

@testset "Termination Conditions: Allocations" begin
    @testset "Mode: $(tcond)" for tcond in TERMINATION_CONDITIONS
        @test_nowarn DiffEqBase.check_convergence(tcond, du, u, uprev, 1e-3, 1e-3)
    end
end
