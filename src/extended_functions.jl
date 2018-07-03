# Implements passing in Jacobians (and possibly other functions) via
# function overloading:
#
# - f(...) - objective function
# - f(Val{:jac}, ...) - Jacobian of objective function
# - the details of what `...` needs to be depends on the
#   AbstractODEProblem subtype

# Method_exists does not work:
#
# julia> f(::Val{:jac}, a, b, c) = 5
# f (generic function with 1 method)
#
# julia> hasmethod(f, Tuple{Val{:jac}, Vararg})
# false
#
# Thus hand-code it:
check_first_arg(f,T::Type) = check_first_arg(typeof(f),T)
function check_first_arg(::Type{F}, T::Type) where F
    typ = Tuple{Any, T, Vararg}
    typ2 = Tuple{Any, Type{T}, Vararg} # This one is required for overloaded types
    method_table = Base.MethodList(F.name.mt) # F.name.mt gets the method-table
    for m in method_table
        (m.sig<:typ || m.sig<:typ2) && return true
    end
    return false
end
# Standard
@traitdef HasJac{F}
@traitdef HastGrad{F}
__has_jac(f) = check_first_arg(f, Val{:jac})
__has_tgrad(f) = check_first_arg(f, Val{:tgrad})
@generated SimpleTraits.trait(::Type{HasJac{F}}) where {F} = __has_jac(F) ? :(HasJac{F}) : :(Not{HasJac{F}})
@generated SimpleTraits.trait(::Type{HastGrad{F}}) where {F} = __has_tgrad(F) ? :(HastGrad{F}) : :(Not{HastGrad{F}})
has_jac(f::T) where {T} = istrait(HasJac{T})
has_tgrad(f::T) where {T} = istrait(HastGrad{T})

# Performance
@traitdef HasExpJac{F}
@traitdef HasInvJac{F}
@traitdef HasInvW{F}
@traitdef HasInvW_t{F}
__has_expjac(f) = check_first_arg(f, Val{:expjac})
__has_invjac(f) = check_first_arg(f, Val{:invjac})
__has_invW(f) = check_first_arg(f, Val{:invW})
__has_invW_t(f) = check_first_arg(f, Val{:invW_t})
@generated SimpleTraits.trait(::Type{HasExpJac{F}}) where {F} = __has_expjac(F) ? :(HasExpJac{F}) : :(Not{HasExpJac{F}})
@generated SimpleTraits.trait(::Type{HasInvJac{F}}) where {F} = __has_invjac(F) ? :(HasInvJac{F}) : :(Not{HasInvJac{F}})
@generated SimpleTraits.trait(::Type{HasInvW{F}}) where {F} = __has_invW(F) ? :(HasInvW{F}) : :(Not{HasInvW{F}})
@generated SimpleTraits.trait(::Type{HasInvW_t{F}}) where {F} = __has_invW_t(F) ? :(HasInvW_t{F}) : :(Not{HasInvW_t{F}})
has_expjac(f::T) where {T} = istrait(HasExpJac{T})
has_invjac(f::T) where {T} = istrait(HasInvJac{T})
has_invW(f::T) where {T} = istrait(HasInvW{T})
has_invW_t(f::T) where {T} = istrait(HasInvW_t{T})

# Hessians
@traitdef HasHes{F}
@traitdef HasInvHes{F}
__has_hes(f) = check_first_arg(f, Val{:hes})
__has_invhes(f) = check_first_arg(f, Val{:invhes})
@generated SimpleTraits.trait(::Type{HasHes{F}}) where {F} = __has_hes(F) ? :(HasHes{F}) : :(Not{HasHes{F}})
@generated SimpleTraits.trait(::Type{HasInvHes{F}}) where {F} = __has_invhes(F) ? :(HasInvHes{F}) : :(Not{HasInvHes{F}})
has_hes(f::T) where {T} = istrait(HasHes{T})
has_invhes(f::T) where {T} = istrait(HasInvHes{T})

# Parameter-Based
@traitdef HasParamDeriv{F}
@traitdef HasParamJac{F}
__has_paramderiv(f) = check_first_arg(f, Val{:deriv})
__has_paramjac(f) = check_first_arg(f, Val{:paramjac})
@generated SimpleTraits.trait(::Type{HasParamDeriv{F}}) where {F} = __has_paramderiv(F) ? :(HasParamDeriv{F}) : :(Not{HasParamDeriv{F}})
@generated SimpleTraits.trait(::Type{HasParamJac{F}}) where {F} = __has_paramjac(F) ? :(HasParamJac{F}) : :(Not{HasParamJac{F}})
has_paramderiv(f::T) where {T} = istrait(HasParamDeriv{T})
has_paramjac(f::T) where {T} = istrait(HasParamJac{T})

## Parameter Names Check
Base.@pure has_syms(f::T) where {T} = :syms ∈ fieldnames(typeof(f))

## Analytical Solution Check
@traitdef HasAnalytic{F}
__has_analytic(f) = check_first_arg(f, Val{:analytic})
@generated SimpleTraits.trait(::Type{HasAnalytic{F}}) where {F} = __has_analytic(F) ? :(HasAnalytic{F}) : :(Not{HasAnalytic{F}})
has_analytic(f::T) where {T} = istrait(HasAnalytic{T})

# now a trait methods can dispatch on this:
# @traitfn fn(g::::HasJac, ...) = ...
# @traitfn fn(g::::(!HasJac), ...) = ...
