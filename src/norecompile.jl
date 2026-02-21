struct OrdinaryDiffEqTag end

const NORECOMPILE_ARGUMENT_MESSAGE = """
No-recompile mode is only supported for state arguments
of type `Vector{Float64}`, time arguments of `Float64`
and parameter arguments of type `Vector{Float64}` or
`SciMLBase.NullParameters`.
"""

struct NoRecompileArgumentError <: Exception
    args::Any
end

function Base.showerror(io::IO, e::NoRecompileArgumentError)
    println(io, NORECOMPILE_ARGUMENT_MESSAGE)
    print(io, "Attempted arguments: ")
    return print(io, e.args)
end

function unwrap_fw(fw::FunctionWrapper)
    return fw.obj[]
end

"""
    DEIIPFunctionWrapper{duType, uType, pType, tType}

Wrapper struct around a `FunctionWrappersWrapper` containing a single `FunctionWrapper`
for an in-place function `f!(du, u, p, t) -> Nothing`.

Compared to a raw `FunctionWrappersWrapper`, this struct exposes only 4 type
parameters (du, u, p, t types) instead of the full nested
`FunctionWrappersWrapper{Tuple{FunctionWrapper{...}}, false}` type, which
significantly reduces type string length in stack traces.

Used by solvers that do **not** require ForwardDiff internally (e.g. Tsit5, Verner).
See also [`DEIIPFunctionWrapperForwardDiff`](@ref) for the ForwardDiff-aware variant.
"""
struct DEIIPFunctionWrapper{duType, uType, pType, tType}
    fw::FunctionWrappersWrappers.FunctionWrappersWrapper{
        Tuple{FunctionWrapper{Nothing, Tuple{duType, uType, pType, tType}}},
        false,
    }
end

(f::DEIIPFunctionWrapper)(args...) = f.fw(args...)
SciMLBase.isfunctionwrapper(::DEIIPFunctionWrapper) = true

"""
    DEIIPFunctionWrapperVF64{pType}

VF64-specialized alias: `DEIIPFunctionWrapper{Vector{Float64}, Vector{Float64}, pType, Float64}`.
Matches the wrapper produced for the common in-place `Vector{Float64}` ODE case
when ForwardDiff is **not** used by the solver.
"""
const DEIIPFunctionWrapperVF64{pType} =
    DEIIPFunctionWrapper{Vector{Float64}, Vector{Float64}, pType, Float64}

"""
    DEIIPFunctionWrapperForwardDiff{T1, T2, T3, T4, dT1, dT2, dT4}

Wrapper struct around a `FunctionWrappersWrapper` containing 4 `FunctionWrapper`
entries for an in-place function `f!(du, u, p, t) -> Nothing` with ForwardDiff support.

The 4 wrappers cover:
1. Base types: `(T1, T2, T3, T4)`
2. Dual state: `(dT1, dT2, T3, T4)`
3. Dual time: `(dT1, T2, T3, dT4)`
4. Dual state+time: `(dT1, dT2, T3, dT4)`

Compared to a raw `FunctionWrappersWrapper`, this struct exposes 7 type parameters
instead of repeating the full `FunctionWrapper{Nothing, Tuple{...}}` 4 times,
significantly reducing type string length in stack traces.

Used by solvers that require ForwardDiff internally (e.g. Rosenbrock, implicit methods).
See also [`DEIIPFunctionWrapper`](@ref) for the simpler non-ForwardDiff variant.
"""
struct DEIIPFunctionWrapperForwardDiff{T1, T2, T3, T4, dT1, dT2, dT4}
    fw::FunctionWrappersWrappers.FunctionWrappersWrapper{
        Tuple{
            FunctionWrapper{Nothing, Tuple{T1, T2, T3, T4}},
            FunctionWrapper{Nothing, Tuple{dT1, dT2, T3, T4}},
            FunctionWrapper{Nothing, Tuple{dT1, T2, T3, dT4}},
            FunctionWrapper{Nothing, Tuple{dT1, dT2, T3, dT4}},
        },
        false,
    }
end

(f::DEIIPFunctionWrapperForwardDiff)(args...) = f.fw(args...)
SciMLBase.isfunctionwrapper(::DEIIPFunctionWrapperForwardDiff) = true

# Union for isa checks (avoid double-wrapping)
const AnyFunctionWrapper = Union{
    FunctionWrappersWrappers.FunctionWrappersWrapper,
    DEIIPFunctionWrapper,
    DEIIPFunctionWrapperForwardDiff,
}

"""
    wrapfun_iip_simple(ff, du, u, p, t)

Wrap an in-place function `ff(du, u, p, t) -> Nothing` into a [`DEIIPFunctionWrapper`](@ref)
(single `FunctionWrapper`, no ForwardDiff support).

Unlike [`wrapfun_iip`](@ref), this function is **not** overridden by the ForwardDiff
extension, so it avoids creating method-table backedges that cause invalidation.
Use this when the solver does not need ForwardDiff internally.
"""
function wrapfun_iip_simple(ff, du, u, p, t)
    inner = FunctionWrappersWrappers.FunctionWrappersWrapper(
        Void(ff), (typeof((du, u, p, t)),), (Nothing,)
    )
    return DEIIPFunctionWrapper(inner)
end

# Default dispatch assumes no ForwardDiff, gets added in the new dispatch
function wrapfun_iip(ff, inputs)
    return FunctionWrappersWrappers.FunctionWrappersWrapper(
        Void(ff), (typeof(inputs),), (Nothing,)
    )
end

function wrapfun_oop(ff, inputs)
    return FunctionWrappersWrappers.FunctionWrappersWrapper(
        ff, (typeof(inputs),), (typeof(inputs[1]),)
    )
end

# Wrap an in-place Jacobian function jac!(J, u, p, t) -> Nothing.
# Unlike the RHS, the Jacobian is not called with Dual numbers
# (the analytical Jacobian IS the derivative), so we only need
# a single FunctionWrapper variant.
function wrapfun_jac_iip(jac_f, inputs)
    return FunctionWrappersWrappers.FunctionWrappersWrapper(
        Void(jac_f), (typeof(inputs),), (Nothing,)
    )
end
