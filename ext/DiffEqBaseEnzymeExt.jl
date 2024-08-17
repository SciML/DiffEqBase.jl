module DiffEqBaseEnzymeExt

using DiffEqBase
import DiffEqBase: value, fastpow
using Enzyme
import Enzyme: Const
using ChainRulesCore

function Enzyme.EnzymeRules.augmented_primal(config::Enzyme.EnzymeRules.ConfigWidth{1},
        func::Const{typeof(DiffEqBase.solve_up)}, ::Type{Duplicated{RT}}, prob,
        sensealg::Union{Const{Nothing}, Const{<:DiffEqBase.AbstractSensitivityAlgorithm}},
        u0, p, args...; kwargs...) where {RT}
    @inline function copy_or_reuse(val, idx)
        if Enzyme.EnzymeRules.overwritten(config)[idx] && ismutable(val)
            return deepcopy(val)
        else
            return val
        end
    end

    @inline function arg_copy(i)
        copy_or_reuse(args[i].val, i + 5)
    end

    res = DiffEqBase._solve_adjoint(
        copy_or_reuse(prob.val, 2), copy_or_reuse(sensealg.val, 3),
        copy_or_reuse(u0.val, 4), copy_or_reuse(p.val, 5),
        SciMLBase.EnzymeOriginator(), ntuple(arg_copy, Val(length(args)))...;
        kwargs...)

    dres = Enzyme.make_zero(res[1])::RT
    tup = (dres, res[2])
    return Enzyme.EnzymeRules.AugmentedReturn{RT, RT, Any}(res[1], dres, tup::Any)
end

function Enzyme.EnzymeRules.reverse(config::Enzyme.EnzymeRules.ConfigWidth{1},
        func::Const{typeof(DiffEqBase.solve_up)}, ::Type{Duplicated{RT}}, tape, prob,
        sensealg::Union{Const{Nothing}, Const{<:DiffEqBase.AbstractSensitivityAlgorithm}},
        u0, p, args...; kwargs...) where {RT}
    dres, clos = tape
    dres = dres::RT
    dargs = clos(dres)
    for (darg, ptr) in zip(dargs, (func, prob, sensealg, u0, p, args...))
        if ptr isa Enzyme.Const
            continue
        end
        if darg == ChainRulesCore.NoTangent()
            continue
        end
        ptr.dval .+= darg
    end
    Enzyme.make_zero!(dres.u)
    return ntuple(_ -> nothing, Val(length(args) + 4))
end

function EnzymeRules.forward(func::Const{typeof(fastpow)},
            RT::Type{<:Union{Const,DuplicatedNoNeed,Duplicated,
            BatchDuplicated,BatchDuplicatedNoNeed}},
            _x::Annotation, _y::Annotation)
    x = _x.val
    y = _y.val
    ret = func.val(x.val, y.val)
    dxval = x.dval * y * (fastpow(x,y - 1))
    dyval = x isa Real && x<=0 ? Base.oftype(float(x), NaN) : y.dval*(fastpow(x,y))*log(x)
    return Duplicated(ret, dxval + dyval)
end

function EnzymeRules.augmented_primal(config::ConfigWidth{1}, 
                          func::Const{typeof(fastpow)}, 
                          ::Type{<:Active},
                          x::Active, x::Active)
    if EnzymeRules.needs_primal(config)
        primal = func.val(x.val, y.val)
    else
        primal = nothing
    end
    return EnzymeRules.AugmentedReturn(primal, nothing, nothing)
end

function EnzymeRules.reverse(config::EnzymeRules.ConfigWidth{1}, 
                            func::Const{DiffEqBase.fastpow}, dret, tape::Nothing,
                            _x, _y)
    x = _x.val
    y = _y.val
    dxval = x.dval * y * (fastpow(x,y - 1))
    dyval = x isa Real && x<=0 ? Base.oftype(float(x), NaN) : y.dval*(fastpow(x,y))*log(x)
    return (dxval, dyval)
end

end
