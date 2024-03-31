using BenchmarkTools, CUDA, DiffEqBase, Test, LinearAlgebra
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
        for nfn in (Base.Fix1(maximum, abs), Base.Fix2(norm, 2), Base.Fix2(norm, Inf))
            tcond = DiffEqBase.set_termination_mode_internalnorm(tcond, nfn)
            @test_nowarn DiffEqBase.check_convergence(tcond, du, u, uprev, 1e-3, 1e-3)
        end
    end
end
