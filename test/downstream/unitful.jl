using Unitful, OrdinaryDiffEq, Test
f(du, u, p, t) = du .= 3 * u"1/s" * u
prob = ODEProblem(f, [2.0u"m"], (0.0u"s", Inf * u"s"))
intg = init(prob, Tsit5())
@test_nowarn step!(intg, 0.02u"s", true)
