using DiffEqBase, Base.Test

# Basic expressions
@test macroexpand(:(@muladd a*b+c)) == :($(Base.muladd)(a, b, c))
@test macroexpand(:(@muladd c+a*b)) == :($(Base.muladd)(a, b, c))
@test macroexpand(:(@muladd b*a+c)) == :($(Base.muladd)(b, a, c))
@test macroexpand(:(@muladd c+b*a)) == :($(Base.muladd)(b, a, c))

# Multiple multiplications
@test macroexpand(:(@muladd a*b+c*d)) == :($(Base.muladd)(a, b, c*d))
@test macroexpand(:(@muladd a*b+c*d+e*f)) == :($(Base.muladd)(a, b,
                                                              $(Base.muladd)(c, d, e*f)))
@test macroexpand(:(@muladd a*(b*c+d)+e)) == :($(Base.muladd)(a,
                                                              $(Base.muladd)(b, c, d), e))

# Dot calls
@test macroexpand(:(@. @muladd a*b+c)) == :($(Base.muladd).(a, b, c))
@test macroexpand(:(@muladd @. a*b+c)) == :($(Base.muladd).(a, b, c))
@test macroexpand(:(@muladd a.*b+c)) == :(a.*b+c)
@test macroexpand(:(@muladd a*b.+c)) == :(a*b.+c)
@test macroexpand(:(@muladd f.(a)*b+c)) == :($(Base.muladd)(f.(a), b, c))
@test macroexpand(:(@muladd a*f.(b)+c)) == :($(Base.muladd)(a, f.(b), c))
@test macroexpand(:(@muladd a*b+f.(c))) == :($(Base.muladd)(a, b, f.(c)))

# Nested expressions
@test macroexpand(:(@muladd f(x, y, z) = x*y+z)) == :(f(x, y, z) = $(Base.muladd)(x, y, z))
@test macroexpand(:(@muladd function f(x, y, z) x*y+z end)) ==
    :(function f(x, y, z) $(Base.muladd)(x, y, z) end)
@test macroexpand(:(@muladd for i in 1:n z = x*i + y end)) ==
    :(for i in 1:n z = $(Base.muladd)(x, i, y) end)
