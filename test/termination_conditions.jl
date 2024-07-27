using BenchmarkTools, DiffEqBase, LinearAlgebra, Test

du = rand(4)
u = rand(4)
uprev = rand(4)

const TERMINATION_CONDITIONS = [
    RelTerminationMode(), NormTerminationMode(), RelNormTerminationMode(),
    AbsTerminationMode(), AbsNormTerminationMode(), RelSafeTerminationMode(),
    AbsSafeTerminationMode(), RelSafeBestTerminationMode(), AbsSafeBestTerminationMode()
]

@testset "Termination Conditions: Allocations" begin
    @testset "Mode: $(tcond)" for tcond in TERMINATION_CONDITIONS
        for nfn in (Base.Fix1(maximum, abs), Base.Fix2(norm, 2), Base.Fix2(norm, Inf))
            tcond = DiffEqBase.set_termination_mode_internalnorm(tcond, nfn)
            @test (@ballocated DiffEqBase.check_convergence($tcond, $du, $u, $uprev, 1e-3,
                1e-3)) == 0
        end
    end
end
