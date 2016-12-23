# Checks that promotion to more general DEProblem types works.
using DiffEqProblemLibrary, Sundials

# prob = prob_ode_2Dlinear this one does not work with IDA
for prob in [prob_ode_linear, prob_ode_vanderpol]

    @test isa(convert(DAEProblem, prob), DAEProblem)
    if isa(prob,ODETestProblem)
        @test isa(promote(prob, IDA()), DAETestProblem)
        @test isa(convert(DAETestProblem, prob), DAETestProblem)
    else
        @test isa(promote(prob, IDA()), DAEProblem)
    end

    sol1 = solve(prob,IDA())
end
