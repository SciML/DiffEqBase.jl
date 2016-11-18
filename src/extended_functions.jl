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
# julia> method_exists(f, Tuple{Val{:jac}, Vararg})
# false
#
# Thus hand-code it:
check_first_arg(f,T::Type) = check_first_arg(typeof(f),T)
function check_first_arg{F}(::Type{F}, T::Type)
    typ = Tuple{Any, T, Vararg}
    for m in Base.MethodList(F.name.mt) # F.name.mt gets the method-table
        m.sig<:typ && return true
    end
    return false
end

# Standard
@traitdef HasJac{F}
@traitdef HastGrad{F}
has_jac(f) = check_first_arg(f, Val{:jac})
has_tgrad(f) = check_first_arg(f, Val{:tgrad})
@generated SimpleTraits.trait{F}(::Type{HasJac{F}}) = has_jac(F) ? :(HasJac{F}) : :(Not{HasJac{F}})
@generated SimpleTraits.trait{F}(::Type{HastGrad{F}}) = has_tgrad(F) ? :(HastGrad{F}) : :(Not{HastGrad{F}})

# Performance
@traitdef HasInvJac{F}
@traitdef HasInvW{F}
@traitdef HasInvW_t{F}
has_invjac(f) = check_first_arg(f, Val{:invjac})
has_invW(f) = check_first_arg(f, Val{:invW})
has_invW_t(f) = check_first_arg(f, Val{:invW_t})
@generated SimpleTraits.trait{F}(::Type{HasInvJac{F}}) = has_invjac(F) ? :(HasInvJac{F}) : :(Not{HasInvJac{F}})
@generated SimpleTraits.trait{F}(::Type{HasInvW{F}}) = has_invW(F) ? :(HasInvW{F}) : :(Not{HasInvW{F}})
@generated SimpleTraits.trait{F}(::Type{HasInvW_t{F}}) = has_invW_t(F) ? :(HasInvW_t{F}) : :(Not{HasInvW_t{F}})

# Hessians
@traitdef HasHes{F}
@traitdef HasInvHes{F}
has_hes(f) = check_first_arg(f, Val{:hes})
has_invhes(f) = check_first_arg(f, Val{:invhes})
@generated SimpleTraits.trait{F}(::Type{HasHes{F}}) = has_hes(F) ? :(HasHes{F}) : :(Not{HasHes{F}})
@generated SimpleTraits.trait{F}(::Type{HasInvHes{F}}) = has_invhes(F) ? :(HasInvHes{F}) : :(Not{HasInvHes{F}})

# Parameter-Based
@traitdef HasParamDeriv{F}
@traitdef HasParamJac{F}
has_paramderiv(f) = check_first_arg(f, Val{:deriv})
has_paramjac(f) = check_first_arg(f, Val{:paramjac})
@generated SimpleTraits.trait{F}(::Type{HasParamDeriv{F}}) = has_paramderiv(F) ? :(HasParamDeriv{F}) : :(Not{HasParamDeriv{F}})
@generated SimpleTraits.trait{F}(::Type{HasParamJac{F}}) = has_paramjac(F) ? :(HasParamJac{F}) : :(Not{HasParamJac{F}})

# now a trait methods can dispatch on this:
# @traitfn fn(g::::HasJac, ...) = ...
# @traitfn fn(g::::(!HasJac), ...) = ...
