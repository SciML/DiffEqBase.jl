module DiffEqBaseChainRulesCoreExt

using DiffEqBase
using DiffEqBase.SciMLBase
import DiffEqBase: numargs, AbstractSensitivityAlgorithm, AbstractDEProblem

import ChainRulesCore
import ChainRulesCore: NoTangent

ChainRulesCore.rrule(::typeof(numargs), f) = (numargs(f), df -> (NoTangent(), NoTangent()))
ChainRulesCore.@non_differentiable DiffEqBase.checkkwargs(kwargshandle)

function ChainRulesCore.frule(
        ::typeof(DiffEqBase.solve_up), prob,
        sensealg::Union{Nothing, AbstractSensitivityAlgorithm},
        u0, p, args...; originator = SciMLBase.ChainRulesOriginator(),
        kwargs...
    )
    return DiffEqBase._solve_forward(
        prob, sensealg, u0, p,
        originator, args...;
        kwargs...
    )
end

function ChainRulesCore.rrule(
        ::typeof(DiffEqBase.solve_up), prob::AbstractDEProblem,
        sensealg::Union{Nothing, AbstractSensitivityAlgorithm},
        u0, p, args...; originator = SciMLBase.ChainRulesOriginator(),
        kwargs...
    )
    return DiffEqBase._solve_adjoint(
        prob, sensealg, u0, p,
        originator, args...;
        kwargs...
    )
end

end
