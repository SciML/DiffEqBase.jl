"""
Defines a differential algebraic equation (DAE) (or implicitly defined
ODE if there are no purely algebraic equations) of form: `f(t,u,du/dt) = 0`

Input:
- `f(t,u,du,out)` or `f(t,u,du)`
- `u0` initial conditions (IC) on `u`
- `du0` IC on `du`
- `tspan` tuple of start and end-time `(t0,tend)`

Extra function, such as the Jacobian can be provided by function overloading:
- `f(::Val{:jac}, t, u, du, alpha, out)` where `jac = df/du + alpha* df/d(du)`

Notes:
- the types of inputs determines the types used in the numerical
  solver.  For instance, using `Int`s in `tspan` will not work.
- ideally the IC are consistent, i.e. `f(tspan[1], u0, du0)==0`.  Some
  solvers will solve for consistent IC but this may fail.
"""
type DAEProblem{uType,duType,tType,isinplace,F} <: AbstractDAEProblem{uType,duType,tType,isinplace,F}
  f::F
  u0::uType
  du0::duType
  tspan::Tuple{tType,tType}
end

type DAETestProblem{uType,duType,tType,isinplace,F} <: AbstractDAETestProblem{uType,duType,tType,isinplace,F}
  f::F
  u0::uType
  du0::duType
  analytic::Base.Callable
  tspan::Tuple{tType,tType}
end

function DAEProblem(f,u0,du0,tspan)
  isinplace = numparameters(f)>=4
  DAEProblem{typeof(u0),typeof(du0),eltype(tspan),isinplace,typeof(f)}(f,u0,du0,tspan)
end

function DAETestProblem(f,u0,du0,analytic,tspan=(0.0,1.0))
  isinplace = numparameters(f)>=4
  DAETestProblem{typeof(u0),typeof(du0),eltype(tspan),isinplace,typeof(f)}(f,u0,du0,analytic,tspan)
end

function print{uType,duType,tType,isinplace,F}(io::IO, prob::AbstractDAEProblem{uType,duType,tType,isinplace,F})
  println(io,"AbstractDAEProblem")
  println(io,"Independent Variable Types: $uType")
  println(io,"Depdendent Variable Type: $tType")
  println(io,"Function is in-place? $isinplace")
  nothing
end

#=
function show{uType,tType,isinplace,F}(io::IO,prob::AbstractODEProblem{uType,tType,isinplace,F})
  println(io,"AbstractDAEProblem{$uType,$tType,$isinplace}")
  nothing
end
=#

###
# Promotion & conversion
promote(prob::ODEProblem, ::AbstractDAEAlgorithm) = convert(DAEProblem, prob)
promote(prob::ODETestProblem, ::AbstractDAEAlgorithm) = convert(DAETestProblem, prob)
function convert{P<:AbstractDAEProblem}(::Type{P}, prob::AbstractODEProblem)
    @unpack f,u0,tspan = prob
    f_new = convert_obj_fn(P,prob)
    # IC
    if isinplace(prob)
        du0 = similar(u0) # TODO: this needs to be Unit-stable (Chris, help!)
        f(prob.tspan[1], u0, du0)
    else
        du0 = f(prob.tspan[1], u0)
    end
    # overloaded functions
    for m in methods_overloaded(f)
        val = get_first_arg(m) # get Val{:jac} etc
        convert_overloaded!(DAEProblem, prob, val, f_new) # updates f_new in-place
    end
    # construct
    if P==DAEProblem
        P(f_new, u0, du0, tspan)
    elseif isa(prob,ODETestProblem) && P==DAETestProblem
        P(f_new, u0, du0, prob.analytic, tspan)
    else
        error("Cannot convert non-test problem into test-problem.")
    end
end

#"Convert the objective function `f` or prob"
@traitfn function convert_obj_fn{T<:AbstractDAEProblem, P2<:AbstractODEProblem}(::Type{T}, prob::P2::IsInplace)
    f = prob.f
    return function (t,u,du,out)
        f(t,u,out)
        out .= .-(out, du) # change in 0.6
    end
end
@traitfn function convert_obj_fn{T<:AbstractDAEProblem, P2<:AbstractODEProblem}(::Type{T}, prob::P2::(!IsInplace))
    f = prob.f
    return (t,u,du) -> .-(f(t,u), du) # change in 0.6
end

"""
Converts one method of a overloaded generic function to the new problem-type.

`convert_overloaded!(T, prob, Val{:sometag}, f_new)`

Input:
- T: type of new AbstractIVPProblem
- prob: old problem instance
- Val{:sometag}: convert `:sometag` method
- `f_new`: the generic function to which the converted method is added

Out:
- nothing
"""
function convert_overloaded!{T<:AbstractIVPProblem}(::Type{T}, prob::AbstractIVPProblem, ::Val, f_new)
    warn("No conversion from $(typeof(prob)) to $T for overloaded function method $val.  Not converting...")
end
@traitfn function convert_overloaded!{T<:AbstractDAEProblem, P2<:AbstractODEProblem,Fnew}(::Type{T}, prob::P2::IsInplace, ::Val{:jac}, f_new::Fnew)
    f = prob.f
    function (::Fnew)(::Val{:jac}, t, u, du, alpha, out)
        f(Val{:jac}(), t, u, out)
        out .= .+(out, alpha*I) # alpha*df/ddu == diag(alpha)
    end
    nothing
end
@traitfn function convert_overloaded!{T<:AbstractDAEProblem, P2<:AbstractODEProblem,Fnew}(::Type{T}, prob::P2::(!IsInplace), ::Val{:jac}, f_new::Fnew)
    f = prob.f
    function (::Fnew)(::Val{:jac}, t, u, du, alpha)
        f(Val{:jac}(), t, u, out) + alpha*I
    end
    nothing
end
