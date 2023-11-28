using BenchmarkTools, DiffEqBase, Test

du = rand(4)
u = rand(4)
uprev = rand(4)

const TERMINATION_CONDITIONS = [
    SteadyStateDiffEqTerminationMode(), SimpleNonlinearSolveTerminationMode(),
    NormTerminationMode(), RelTerminationMode(), RelNormTerminationMode(),
    AbsTerminationMode(), AbsNormTerminationMode(), RelSafeTerminationMode(),
    AbsSafeTerminationMode(), RelSafeBestTerminationMode(), AbsSafeBestTerminationMode(),
]

@testset "Termination Conditions: Allocations" begin
    @testset "Mode: $(tcond)" for tcond in TERMINATION_CONDITIONS
        @test (@ballocated DiffEqBase.check_convergence($tcond, $du, $u, $uprev, 1e-3,
            1e-3)) == 0 
    end
end
