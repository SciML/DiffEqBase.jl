module DiffEqBaseEnzymeExt

@static if isempty(VERSION.prerelease)
    using DiffEqBase
    import DiffEqBase: value
    using Enzyme
    import Enzyme: Const
    using ChainRulesCore

    function Enzyme.EnzymeRules.augmented_primal(
            config::Enzyme.EnzymeRules.RevConfigWidth{1},
            func::Const{typeof(DiffEqBase.solve_up)}, ::Type{RT}, prob,
            sensealg::Union{
                Const{Nothing}, Const{<:DiffEqBase.AbstractSensitivityAlgorithm}},
            u0, p, args...; kwargs...) where {RT <: Union{Duplicated, DuplicatedNoNeed}}
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

        ResType = typeof(res[1])
        dres = Enzyme.make_zero(res[1])::ResType
        tup = (dres, res[2])
        return Enzyme.EnzymeRules.AugmentedReturn{ResType, ResType, Any}(res[1], dres, tup::Any)
    end

    function Enzyme.EnzymeRules.reverse(config::Enzyme.EnzymeRules.RevConfigWidth{1},
            func::Const{typeof(DiffEqBase.solve_up)}, ::Type{RT}, tape, prob,
            sensealg::Union{
                Const{Nothing}, Const{<:DiffEqBase.AbstractSensitivityAlgorithm}},
            u0, p, args...; kwargs...) where {RT <: Union{Duplicated, DuplicatedNoNeed}}
        dres, clos = tape
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
end

end
