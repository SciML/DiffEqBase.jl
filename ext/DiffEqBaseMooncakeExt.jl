module DiffEqBaseMooncakeExt

using DiffEqBase, Mooncake
using DiffEqBase: SciMLBase
using SciMLBase: ADOriginator, MooncakeOriginator
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

# Dispatch for auto-alg
Mooncake.@from_rrule(
    Mooncake.MinimalCtx,
    Tuple{
        typeof(DiffEqBase.solve_up),
        DiffEqBase.AbstractDEProblem,
        Union{Nothing, DiffEqBase.AbstractSensitivityAlgorithm},
        Any,
        Any,
    },
    true,
    )

Mooncake.@zero_adjoint Mooncake.MinimalCtx Tuple{typeof(DiffEqBase.numargs), Any}
Mooncake.@mooncake_overlay DiffEqBase.set_mooncakeoriginator_if_mooncake(x::ADOriginator) = MooncakeOriginator

end