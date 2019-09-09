struct LinearProblem{F} <: DEProblem
    A::F
    b
    u0
end

struct NonlinearProblem{F} <: DEProblem
    f::F
    u0
end

struct QuadratureProblem <: DEProblem
    f::F
    lb
    ub
end
