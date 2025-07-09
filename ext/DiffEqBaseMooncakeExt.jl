module DiffEqBaseMooncakeExt

using DiffEqBase, Mooncake
using DiffEqBase: SciMLBase
using SciMLBase: ADOriginator, ChainRulesOriginator, MooncakeOriginator
import Mooncake: rrule!!, CoDual, zero_fcodual, @is_primitive,
                 @from_rrule, @zero_adjoint, @mooncake_overlay, MinimalCtx,
                 NoPullback

@from_rrule(MinimalCtx,
    Tuple{
        typeof(DiffEqBase.solve_up),
        DiffEqBase.AbstractDEProblem,
        Union{Nothing, DiffEqBase.AbstractSensitivityAlgorithm},
        Any,
        Any,
        Any
    },
    true,)

# Dispatch for auto-alg
@from_rrule(MinimalCtx,
    Tuple{
        typeof(DiffEqBase.solve_up),
        DiffEqBase.AbstractDEProblem,
        Union{Nothing, DiffEqBase.AbstractSensitivityAlgorithm},
        Any,
        Any
    },
    true,)

@zero_adjoint MinimalCtx Tuple{typeof(DiffEqBase.numargs), Any}
@is_primitive MinimalCtx Tuple{
    typeof(DiffEqBase.set_mooncakeoriginator_if_mooncake), SciMLBase.ChainRulesOriginator
}

@mooncake_overlay DiffEqBase.set_mooncakeoriginator_if_mooncake(x::SciMLBase.ADOriginator) = SciMLBase.MooncakeOriginator()

function rrule!!(
        f::CoDual{typeof(DiffEqBase.set_mooncakeoriginator_if_mooncake)},
        X::CoDual{SciMLBase.ChainRulesOriginator}
)
    return zero_fcodual(SciMLBase.MooncakeOriginator()), NoPullback(f, X)
end

end
