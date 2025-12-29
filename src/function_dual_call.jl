# Dual in-place/out-of-place call support for ODE functions
# Allows calling in-place functions with out-of-place syntax and vice versa
# See: https://github.com/SciML/DiffEqBase.jl/issues/361

using SciMLBase: ODEFunction, SplitFunction, DiscreteFunction, DAEFunction,
                 DDEFunction, SDEFunction, SDDEFunction, ImplicitDiscreteFunction,
                 NonlinearFunction, AbstractSciMLOperator, numargs

# Helper to check if wrapped function is an operator (which has different call conventions)
_is_operator(f) = f.f isa AbstractSciMLOperator

# Helper to get the number of arguments the wrapped function expects
# Returns the first element of numargs result (handling multiple methods)
function _numargs(f)
    nargs = numargs(f)
    return nargs isa AbstractVector ? first(nargs) : nargs
end

# =============================================================================
# ODEFunction: standard signature is (u, p, t) for OOP and (du, u, p, t) for IIP
#
# NOTE: DynamicalODEFunction stores internal ODEFunctions with non-standard
# signatures: (v, u, p, t) for OOP and (dv, v, u, p, t) for IIP.
# We use numargs to detect and avoid interfering with these.
#
# Standard ODEFunction:
#   - OOP: f.f has 3 args (u, p, t)
#   - IIP: f.f has 4 args (du, u, p, t)
#
# DynamicalODEFunction internal:
#   - OOP: f.f has 4 args (v, u, p, t)
#   - IIP: f.f has 5 args (dv, v, u, p, t)
# =============================================================================

# IIP function called with OOP style: (u, p, t) -> allocate du and return it
# Only applies to standard ODEFunction (4-arg wrapped function), not Dynamical internal (5-arg)
function (f::ODEFunction{true})(u, p, t)
    if _is_operator(f)
        # Operators have different call conventions, let them handle it
        return f.f(u, u, p, t)
    end
    # Check that this is a standard IIP function (4 args), not a DynamicalODE internal (5 args)
    if _numargs(f.f) != 4
        # Fall through to the default behavior in SciMLBase
        return f.f(u, p, t)
    end
    du = similar(u)
    f.f(du, u, p, t)
    return du
end

# OOP function called with IIP style: (du, u, p, t) -> fill du with result
# Only applies to standard ODEFunction (3-arg wrapped function), not Dynamical internal (4-arg)
function (f::ODEFunction{false})(du, u, p, t)
    if _is_operator(f)
        # Operators have different call conventions, let them handle it
        f.f(du, u, u, p, t)
        return nothing
    end
    # Check that this is a standard OOP function (3 args), not a DynamicalODE internal (4 args)
    if _numargs(f.f) != 3
        # This is a DynamicalODE internal function, call it with its native 4-arg signature
        return f.f(du, u, p, t)
    end
    du .= f.f(u, p, t)
    return nothing
end

# =============================================================================
# SplitFunction: same signature as ODEFunction
# =============================================================================

# IIP function called with OOP style
function (f::SplitFunction{true})(u, p, t)
    du = similar(u)
    tmp = similar(u)
    f.f1(tmp, u, p, t)
    f.f2(du, u, p, t)
    du .+= tmp
    return du
end

# OOP function called with IIP style
function (f::SplitFunction{false})(du, u, p, t)
    du .= f.f1(u, p, t) + f.f2(u, p, t)
    return nothing
end

# =============================================================================
# DiscreteFunction: signature is (u, p, t) for OOP and (du, u, p, t) for IIP
# =============================================================================

# IIP function called with OOP style
function (f::DiscreteFunction{true})(u, p, t)
    du = similar(u)
    f.f(du, u, p, t)
    return du
end

# OOP function called with IIP style
function (f::DiscreteFunction{false})(du, u, p, t)
    du .= f.f(u, p, t)
    return nothing
end

# =============================================================================
# ImplicitDiscreteFunction: signature is (u, p, t) for OOP and (du, u, p, t) for IIP
# =============================================================================

# IIP function called with OOP style
function (f::ImplicitDiscreteFunction{true})(u, p, t)
    du = similar(u)
    f.f(du, u, p, t)
    return du
end

# OOP function called with IIP style
function (f::ImplicitDiscreteFunction{false})(du, u, p, t)
    du .= f.f(u, p, t)
    return nothing
end

# =============================================================================
# SDEFunction: same signature as ODEFunction
# =============================================================================

# IIP function called with OOP style
function (f::SDEFunction{true})(u, p, t)
    if _is_operator(f)
        return f.f(u, u, p, t)
    end
    # Check that this is a standard IIP function (4 args)
    if _numargs(f.f) != 4
        return f.f(u, p, t)
    end
    du = similar(u)
    f.f(du, u, p, t)
    return du
end

# OOP function called with IIP style
function (f::SDEFunction{false})(du, u, p, t)
    if _is_operator(f)
        f.f(du, u, u, p, t)
        return nothing
    end
    # Check that this is a standard OOP function (3 args)
    if _numargs(f.f) != 3
        return f.f(du, u, p, t)
    end
    du .= f.f(u, p, t)
    return nothing
end

# =============================================================================
# DAEFunction: signature is (du, u, p, t) for OOP and (out, du, u, p, t) for IIP
# =============================================================================

# IIP function called with OOP style
function (f::DAEFunction{true})(du, u, p, t)
    # Standard IIP DAE has 5 args
    if _numargs(f.f) != 5
        return f.f(du, u, p, t)
    end
    out = similar(u)
    f.f(out, du, u, p, t)
    return out
end

# OOP function called with IIP style
function (f::DAEFunction{false})(out, du, u, p, t)
    # Standard OOP DAE has 4 args
    if _numargs(f.f) != 4
        return f.f(out, du, u, p, t)
    end
    out .= f.f(du, u, p, t)
    return nothing
end

# =============================================================================
# DDEFunction: signature is (u, h, p, t) for OOP and (du, u, h, p, t) for IIP
#
# NOTE: DynamicalDDEFunction stores internal DDEFunctions with non-standard
# signatures: (v, u, h, p, t) for OOP and (dv, v, u, h, p, t) for IIP.
#
# Standard DDEFunction:
#   - OOP: f.f has 4 args (u, h, p, t)
#   - IIP: f.f has 5 args (du, u, h, p, t)
#
# DynamicalDDEFunction internal:
#   - OOP: f.f has 5 args (v, u, h, p, t)
#   - IIP: f.f has 6 args (dv, v, u, h, p, t)
# =============================================================================

# IIP function called with OOP style
function (f::DDEFunction{true})(u, h, p, t)
    # Standard IIP DDE has 5 args
    if _numargs(f.f) != 5
        return f.f(u, h, p, t)
    end
    du = similar(u)
    f.f(du, u, h, p, t)
    return du
end

# OOP function called with IIP style
function (f::DDEFunction{false})(du, u, h, p, t)
    # Standard OOP DDE has 4 args
    if _numargs(f.f) != 4
        return f.f(du, u, h, p, t)
    end
    du .= f.f(u, h, p, t)
    return nothing
end

# =============================================================================
# SDDEFunction: same signature as DDEFunction
# =============================================================================

# IIP function called with OOP style
function (f::SDDEFunction{true})(u, h, p, t)
    # Standard IIP SDDE has 5 args
    if _numargs(f.f) != 5
        return f.f(u, h, p, t)
    end
    du = similar(u)
    f.f(du, u, h, p, t)
    return du
end

# OOP function called with IIP style
function (f::SDDEFunction{false})(du, u, h, p, t)
    # Standard OOP SDDE has 4 args
    if _numargs(f.f) != 4
        return f.f(du, u, h, p, t)
    end
    du .= f.f(u, h, p, t)
    return nothing
end

# =============================================================================
# NonlinearFunction: signature is (u, p) for OOP and (du, u, p) for IIP
# =============================================================================

# IIP function called with OOP style
function (f::NonlinearFunction{true})(u, p)
    du = similar(u)
    f.f(du, u, p)
    return du
end

# OOP function called with IIP style
function (f::NonlinearFunction{false})(du, u, p)
    du .= f.f(u, p)
    return nothing
end
