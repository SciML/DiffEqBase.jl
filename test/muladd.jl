using DiffEqBase, Base.Test

# Basic expressions
@test macroexpand(:(@muladd a*b+c)) == Expr(:call, Base.muladd, :a, :b, :c)
@test macroexpand(:(@muladd c+a*b)) == Expr(:call, Base.muladd, :a, :b, :c)
@test macroexpand(:(@muladd b*a+c)) == Expr(:call, Base.muladd, :b, :a, :c)
@test macroexpand(:(@muladd c+b*a)) == Expr(:call, Base.muladd, :b, :a, :c)

# Multiple multiplications
@test macroexpand(:(@muladd a*b+c*d)) == Expr(:call, Base.muladd, :a, :b,
                                              Expr(:call, :*, :c, :d))
@test macroexpand(:(@muladd a*b+c*d+e*f)) == Expr(:call, Base.muladd, :a, :b,
                                                  Expr(:call, Base.muladd, :c, :d,
                                                       Expr(:call, :*, :e, :f)))
@test macroexpand(:(@muladd a*(b*c+d)+e)) == Expr(:call, Base.muladd, :a,
                                                  Expr(:call, Base.muladd, :b, :c, :d), :e)

# Dot calls
@test macroexpand(:(@. @muladd a*b+c)) == Expr(:., Base.muladd, Expr(:tuple, :a, :b, :c))
@test macroexpand(:(@muladd @. a*b+c)) == Expr(:., Base.muladd, Expr(:tuple, :a, :b, :c))
@test macroexpand(:(@muladd a.*b+c)) == Expr(:., Base.muladd, Expr(:tuple, :a, :b, :c))
@test macroexpand(:(@muladd a*b.+c)) == Expr(:., Base.muladd, Expr(:tuple, :a, :b, :c))
@test macroexpand(:(@muladd f.(a)*b+c)) == Expr(:., Base.muladd,
                                                Expr(:tuple,
                                                     Expr(:., :f, Expr(:tuple, :a)),
                                                     :b, :c))
@test macroexpand(:(@muladd a*f.(b)+c)) == Expr(:., Base.muladd,
                                                Expr(:tuple, :a,
                                                     Expr(:., :f, Expr(:tuple, :b)),
                                                     :c))
@test macroexpand(:(@muladd a*b+f.(c))) == Expr(:., Base.muladd,
                                                Expr(:tuple, :a, :b,
                                                     Expr(:., :f, Expr(:tuple, :c))))
