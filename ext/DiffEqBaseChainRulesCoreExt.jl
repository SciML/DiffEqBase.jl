module DiffEqBaseChainRulesCoreExt

using DiffEqBase
import DiffEqBase: numargs

import ChainRulesCore
import ChainRulesCore: NoTangent

ChainRulesCore.rrule(::typeof(numargs), f) = (numargs(f), df -> (NoTangent(), NoTangent()))
ChainRulesCore.@non_differentiable checkkwargs(kwargshandle)

function ChainRulesCore.frule(::typeof(solve_up), prob,
    sensealg::Union{Nothing, AbstractSensitivityAlgorithm},
    u0, p, args...;
    kwargs...)
    _solve_forward(prob, sensealg, u0, p, SciMLBase.ChainRulesOriginator(), args...;
        kwargs...)
end

function ChainRulesCore.rrule(::typeof(solve_up), prob::AbstractDEProblem,
    sensealg::Union{Nothing, AbstractSensitivityAlgorithm},
    u0, p, args...;
    kwargs...)
    _solve_adjoint(prob, sensealg, u0, p, SciMLBase.ChainRulesOriginator(), args...;
        kwargs...)
end

end