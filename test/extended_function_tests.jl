using Test, DiffEqBase

f123(x,y) = 1
f123(::Val{:jac}, x) = 2

g123(x,y) = 1
g123(::Val{:jac}, x::T) where {T} = 2

h123(x,y) = 1
h123(x::T) where {T} = 2

struct T123 end
(::T123)(a) = 1
(::T123)(::Val{:jac}, a) = 1
t123 = T123()

struct G123 end
(::G123)(a) = 1
(::G123)(b, a) = 1

@test DiffEqBase.__has_jac(f123)
@test DiffEqBase.__has_jac(g123)
@test !DiffEqBase.__has_jac(h123)
@test DiffEqBase.__has_jac(T123())
@test !DiffEqBase.__has_jac(G123())
