# Standard
__has_jac(f) = isdefined(f, :jac)
__has_tgrad(f) = isdefined(f, :tgrad)

# Performance
__has_invW(f) = isdefined(f, :invW)
__has_invW_t(f) = isdefined(f, :invW_t)
has_invW(f::T) where {T} = istrait(HasInvW{T})
has_invW_t(f::T) where {T} = istrait(HasInvW_t{T})

# Parameter-Based
__has_paramderiv(f) = isdefined(f, :deriv)
__has_paramjac(f) = isdefined(f, :paramjac)

## Parameter Names Check
__has_syms(f) = isdefined(f, :syms)

## Analytical Solution Check
__has_analytic(f) = isdefined(f, :analytic)

# Color vector for sparse Jacobian
__has_colorvec(f) = isdefined(f,:colorvec)
