# Tools for registering algorithms in a central place

@enum Grade notgreat good verybest
@enum Class stiff nonstiff
@enum Precision low medium high

immutable AlgSpecs{Alg<:AbstractODEAlgorithm}
    name::String # a human-readable name, e.g. "Radau II/A"
    pkg::Module
    grade::Grade
    class::Class
    precision::Precision
end
getalg{Alg}(::AlgSpecs{Alg}) = Alg

"""
This dictionary should hold all
"""
const ode_algs = Dict{Any, AlgSpecs}()

function add_alg!{Alg<:AbstractODEAlgorithm}(::Type{Alg},
                                             name::String,
                                             grade::Grade,
                                             class::Class,
                                             precision::Precision)
    ode_algs[Alg] = AlgSpecs{Alg}(name, Alg.name.module, grade, class, precision)
end

const dummy_alg = AlgSpecs{AbstractODEAlgorithm}("dummy", DiffEqBase, notgreat, nonstiff, low)

"""
Returns the best algorithm in its (class,precision)
"""
function get_best_alg(class::Class, precision::Precision)
    best = dummy_alg
    for (Alg,v) in ode_algs
        # class and precision need to match
        (class==v.class && precision==v.precision) || continue
        if best.grade>=v.grade
            best = v
        end
    end
    best==dummy_alg && error("No best algorithm found")

    return best
end

# Default solver choice:
function solve(prob::AbstractODEProblem; class=error("For automatic solver selection specify: stiff or nonstiff"),
               precision=error("For automatic solver selection specify: precision (low, medium, high)"),
               kwargs...)
    alg = get_best_alg(class, precision)
    println("Selected algorithm: $(alg.name)")
    solve(prob, getalg(alg); kwargs...)
end
