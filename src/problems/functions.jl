# Functions for AbstractODEProblem and possibly other DEProblem types
#
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

has_jac(f) = check_first_arg(f, Val{:jac})

# Using SimpleTraits to dispatch on existence of jac-method
using SimpleTraits
@traitdef HasJac{F}
@generated SimpleTraits.trait{F}(::Type{HasJac{F}}) = has_jac(F) ? :(HasJac{F}) : :(Not{HasJac{F}})
# now a trait methods can dispatch on this:
# @traitfn fn(g::::HasJac, ...) = ...
# @traitfn fn(g::::(!HasJac), ...) = ...


######################
## Tests (TODO move to test/...)
#####################
f123(x,y) = 1
f123(::Val{:jac}, x) = 2

g123(x,y) = 1
g123{T}(::Val{:jac}, x::T) = 2

h123(x,y) = 1
h123{T}(x::T) = 2

immutable T123 end
(::T123)(a) = 1
(::T123)(::Val{:jac}, a) = 1

immutable G123 end
(::G123)(a) = 1
(::G123)(b, a) = 1

@assert has_jac(f123)
@assert has_jac(g123)
@assert !has_jac(h123)
@assert has_jac(T123())
@assert !has_jac(G123())
# Arguably below could be ==false as T123 is a constructor
# function. However, I'm not sure how to implement this.
@assert has_jac(T123)
@assert !has_jac(G123)

@traitfn tfn(f::::HasJac) = true
@traitfn tfn(f::::(!HasJac)) = false

@assert tfn(f123)
@assert tfn(g123)
@assert !tfn(h123)
@assert tfn(T123())
@assert !tfn(G123())
@assert !tfn(T123)  # NOTE this is inconsistent as has_jac(T123)==true
@assert !tfn(G123)

# ParameterizedFunction:
type  LotkaVolterra <: ParameterizedFunction
         a::Float64
         b::Float64
end
(p::LotkaVolterra)(t,u,du) = begin
         du[1] = p.a * u[1] - p.b * u[1]*u[2]
         du[2] = -3 * u[2] + u[1]*u[2]
end
(p::LotkaVolterra)(::Val{:jac}, t, u, J) = 1

lv = LotkaVolterra(0.0,0.0)

@assert has_jac(lv)
@assert tfn(lv)
