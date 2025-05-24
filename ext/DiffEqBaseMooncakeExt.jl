module DiffEqBaseMooncakeExt

using DiffEqBase, Mooncake
using DiffEqBase: SciMLBase
using SciMLBase: ADOriginator, MooncakeOriginator, ChainRulesOriginator
Mooncake.@from_rrule(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(DiffEqBase.solve_up),
        DiffEqBase.AbstractDEProblem,
        Union{Nothing, DiffEqBase.AbstractSensitivityAlgorithm},
        Any,
        Any,
        Any,
    },
    true,
    )

Mooncake.@zero_adjoint Mooncake.MinimalCtx Tuple{typeof(DiffEqBase.numargs), Any}
Mooncake.@mooncake_overlay DiffEqBase.set_mooncakeoriginator_if_mooncake(x::ChainRulesOriginator) = MooncakeOriginator()

end