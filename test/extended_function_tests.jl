using Base.Test, DiffEqBase, SimpleTraits

f123(x,y) = 1
f123(::Val{:jac}, x) = 2
@test istrait(HasJac{typeof(f123)})

# @code_llvm istrait(HasJac{typeof(f123)})

g123(x,y) = 1
g123{T}(::Val{:jac}, x::T) = 2

h123(x,y) = 1
h123{T}(x::T) = 2

struct T123 end
(::T123)(a) = 1
(::T123)(::Val{:jac}, a) = 1
t123 = T123()

struct G123 end
(::G123)(a) = 1
(::G123)(b, a) = 1

@test has_jac(f123)
@test has_jac(g123)
@test !has_jac(h123)
@test has_jac(T123())
@test !has_jac(G123())

@traitfn tfn(f::::HasJac) = true
@traitfn tfn(f::::(!HasJac)) = false

@traitfn testing(f::::HasJac)   =2
@traitfn testing(f::::!(HasJac))=3.0

@inferred testing(f123)
@inferred testing(g123)
@inferred testing(h123)

@test tfn(f123)
@test tfn(g123)
@test !tfn(h123)
@test tfn(T123())
@test !tfn(G123())
@test !tfn(T123)  # NOTE this is inconsistent as has_jac(T123)==true
@test !tfn(G123)


using ParameterizedFunctions

type  LotkaVolterra <: AbstractParameterizedFunction{true}
          a::Float64
          b::Float64
end
(p::LotkaVolterra)(t,u,du) = begin
        du[1] = p.a * u[1] - p.b * u[1]*u[2]
        du[2] = -3 * u[2] + u[1]*u[2]
end
(p::LotkaVolterra)(::Val{:jac}, t, u, J) = 1

lv = LotkaVolterra(0.0,0.0)

@test has_jac(lv)
@test !has_syms(lv)
@test tfn(lv)

lotka = @ode_def_nohes LotkaVolterraTest begin
  dx = a*x - b*x*y
  dy = -c*y + d*x*y
end a=>1.5 b=>1 c=3 d=1

@test has_syms(lotka)
