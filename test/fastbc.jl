using DiffEqBase: @..
using DiffEqBase
using Test
using InteractiveUtils
using MuladdMacro
function foo9(a, b, c, d, e, f, g, h, i)
    a = DiffEqBase.diffeqbc(a)
    @.. a = b + 0.1 * (0.2c + 0.3d + 0.4e + 0.5f + 0.6g + 0.6h + 0.6i)
    nothing
end
function foo26(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z)
    @muladd @.. a = b + 0.1 * (0.2c + 0.3d + 0.4e + 0.5f + 0.6g + 0.6h + 0.6i + 0.6j + 0.6k + 0.6l + 0.6m + 0.6n + 0.6o + 0.6p + 0.6q + 0.6r + 0.6s + 0.6t + 0.6u + 0.6v + 0.6w + 0.6x + 0.6y + 0.6z)
    nothing
end
a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z = [rand(1000) for i in 1:26];
@allocated foo9(a, b, c, d, e, f, g, h, i)
@allocated foo26(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z)
@test @allocated(foo9(a, b, c, d, e, f, g, h, i)) == 0
@test @allocated(foo26(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, v, w, x, y, z)) == 0
io = IOBuffer()
InteractiveUtils.code_llvm(io, foo9, (Base.typesof)(a, b, c, d, e, f, g, h, i))
str = String(take!(io))
@test occursin("vector.body", str)
