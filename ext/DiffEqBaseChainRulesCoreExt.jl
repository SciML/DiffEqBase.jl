module DiffEqBaseChainRulesCoreExt

using DiffEqBase
import DiffEqBase: numargs

import ChainRulesCore
import ChainRulesCore: NoTangent

ChainRulesCore.rrule(::typeof(numargs), f) = (numargs(f), df -> (NoTangent(), NoTangent()))

end