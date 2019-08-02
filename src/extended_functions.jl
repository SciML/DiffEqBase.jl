# Standard
__has_jac(f) = isdefined(f, :jac)
__has_tgrad(f) = isdefined(f, :tgrad)

# Performance
__has_Wfact(f) = isdefined(f, :Wfact)
__has_Wfact_t(f) = isdefined(f, :Wfact_t)

# Parameter-Based
__has_paramderiv(f) = isdefined(f, :deriv)
__has_paramjac(f) = isdefined(f, :paramjac)

## Parameter Names Check
__has_syms(f) = isdefined(f, :syms)

## Analytical Solution Check
__has_analytic(f) = isdefined(f, :analytic)

# Color vector for sparse Jacobian
__has_color(f) = isdefined(f,:color)
