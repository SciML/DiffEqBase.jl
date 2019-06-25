using Test, DiffEqBase
using DiffEqBase: __has_jac, __has_tgrad, __has_invW, __has_invW_t, __has_paramderiv,
  __has_paramjac, __has_syms, __has_analytic, __has_colorvec

struct Foo
  jac
  tgrad
  invW
  invW_t
  deriv
  paramjac
  syms
  analytic
  colorvec
end

f = Foo(1,1,1,1,1,1,1,1,1)

@test __has_jac(f)
@test __has_tgrad(f)
@test __has_invW(f)
@test __has_invW_t(f)
@test __has_paramderiv(f)
@test __has_paramjac(f)
@test __has_syms(f)
@test __has_analytic(f)
@test __has_colorvec(f)
