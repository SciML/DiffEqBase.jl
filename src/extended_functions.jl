# Standard
@deprecate __has_jac(f) isdefined(f, :jac)
@deprecate __has_tgrad(f) isdefined(f, :tgrad)

# Performance
@deprecate __has_Wfact(f) isdefined(f, :Wfact)
@deprecate __has_Wfact_t(f) isdefined(f, :Wfact_t)

# Parameter-Based
@deprecate __has_paramderiv(f) isdefined(f, :deriv)
@deprecate __has_paramjac(f) isdefined(f, :paramjac)

## Parameter Names Check
@deprecate __has_syms(f) isdefined(f, :syms)

## Analytical Solution Check
@deprecate __has_analytic(f) isdefined(f, :analytic)

# Color vector for sparse Jacobian
@deprecate __has_colorvec(f) isdefined(f,:colorvec)
