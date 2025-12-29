using DiffEqBase, Test, RecursiveArrayTools
import SciMLBase: AbstractSciMLFunction

macro iop_def(funcdef::Expr)
    """Define in- and out-of-place functions simultaneously.

    Call on oop function definition, defines two functions with suffixes _op and _ip.
    """
    @assert funcdef.head ∈ (:function, :(=)) && funcdef.args[1].head == :call

    fname = funcdef.args[1].args[1]

    opname = Symbol("$(fname)_op")
    ipname = Symbol("$(fname)_ip")

    opdef = deepcopy(funcdef)
    opdef.args[1].args[1] = opname

    return quote
        $(esc(opdef))
        $(esc(ipname))(du, args...) = du .= $(esc(opname))(args...)
    end
end

function test_inplace(du, expected, f::AbstractSciMLFunction, args...)
    """Test the in-place version of a function."""
    fill!(du, NaN)
    f(du, args...)
    @test du == expected
end

# Allocate du automatically based on type of expected result
function test_inplace(expected, f::AbstractSciMLFunction, args...)
    test_inplace(similar(expected), expected, f, args...)
end

function test_iop(expected, f_op::AbstractSciMLFunction, f_ip::AbstractSciMLFunction, args...)
    """Test in- and out-of-place version of function both match expected value."""
    @test f_op(args...) == expected
    test_inplace(expected, f_ip, args...)
end

@iop_def f(u, p, t) = p[1] .* u
u = [1.0, 2.0, 3.0]
p = [2.0]
t = 0.0

@testset "ODEFunction with default recompile flag" begin
    odefun = ODEFunction{false}(f_op)
    odefun_ip = ODEFunction{true}(f_ip)
    expected = f_op(u, p, t)
    test_iop(expected, odefun, odefun_ip, u, p, t)
end

@testset "ODEFunction with recompile flag: $rflag" for rflag in (true, false)
    odefun = ODEFunction{false, rflag}(f_op)
    odefun_ip = ODEFunction{true, rflag}(f_ip)
    expected = f_op(u, p, t)
    test_iop(expected, odefun, odefun_ip, u, p, t)
end

# SplitFunction
@iop_def f2(u, p, t) = u .^ 2
sfun = SplitFunction{false}(f_op, f2_op)
sfun_ip = SplitFunction{true}(f_ip, f2_ip; _func_cache = similar(u))
expected = f_op(u, p, t) + f2_op(u, p, t)
test_iop(expected, sfun, sfun_ip, u, p, t)

# DynamicalODEFunction
@iop_def dode_f1(v, u, p, t) = -u
@iop_def dode_f2(v, u, p, t) = p[1] .* v
dodefun = DynamicalODEFunction{false}(dode_f1_op, dode_f2_op)
dodefun_ip = DynamicalODEFunction{true}(dode_f1_ip, dode_f2_ip)
v = [4.0, 5.0, 6.0]
expected = ArrayPartition(dode_f1_op(v, u, p, t), dode_f2_op(v, u, p, t))
test_iop(expected, dodefun, dodefun_ip, ArrayPartition(v, u), p, t)

# DiscreteFunction
dfun = DiscreteFunction{false}(f_op)
dfun_ip = DiscreteFunction{true}(f_ip)
test_iop(f_op(u, p, t), dfun, dfun_ip, u, p, t)

# Type stability
f_analytic(u, p, t) = u
jac = (u, p, t) -> 1
@inferred ODEFunction{false}(f_op, jac = jac)
@inferred DiscreteFunction{false}(f_op, analytic = f_analytic)

# =============================================================================
# Test dual in-place/out-of-place call support (Issue #361)
# Allows calling in-place functions with out-of-place syntax and vice versa
# =============================================================================

@testset "Dual in-place/out-of-place call support" begin
    u = [1.0, 2.0, 3.0]
    p = [2.0]
    t = 0.0

    @testset "ODEFunction" begin
        f_oop(u, p, t) = p[1] .* u
        f_iip(du, u, p, t) = (du .= p[1] .* u)

        ode_oop = ODEFunction{false}(f_oop)
        ode_iip = ODEFunction{true}(f_iip)

        expected = p[1] .* u

        # Standard calls
        @test ode_oop(u, p, t) == expected
        du = zeros(3)
        ode_iip(du, u, p, t)
        @test du == expected

        # Cross calls: IIP function with OOP style
        @test ode_iip(u, p, t) == expected

        # Cross calls: OOP function with IIP style
        du = zeros(3)
        ode_oop(du, u, p, t)
        @test du == expected
    end

    @testset "SplitFunction" begin
        f1_oop(u, p, t) = u
        f2_oop(u, p, t) = u .^ 2
        f1_iip(du, u, p, t) = (du .= u)
        f2_iip(du, u, p, t) = (du .= u .^ 2)

        split_oop = SplitFunction{false}(f1_oop, f2_oop)
        split_iip = SplitFunction{true}(f1_iip, f2_iip; _func_cache = similar(u))

        expected = u + u .^ 2

        # Standard calls
        @test split_oop(u, p, t) == expected
        du = zeros(3)
        split_iip(du, u, p, t)
        @test du == expected

        # Cross calls: IIP function with OOP style
        @test split_iip(u, p, t) == expected

        # Cross calls: OOP function with IIP style
        du = zeros(3)
        split_oop(du, u, p, t)
        @test du == expected
    end

    @testset "DiscreteFunction" begin
        f_oop(u, p, t) = p[1] .* u
        f_iip(du, u, p, t) = (du .= p[1] .* u)

        disc_oop = DiscreteFunction{false}(f_oop)
        disc_iip = DiscreteFunction{true}(f_iip)

        expected = p[1] .* u

        # Standard calls
        @test disc_oop(u, p, t) == expected
        du = zeros(3)
        disc_iip(du, u, p, t)
        @test du == expected

        # Cross calls: IIP function with OOP style
        @test disc_iip(u, p, t) == expected

        # Cross calls: OOP function with IIP style
        du = zeros(3)
        disc_oop(du, u, p, t)
        @test du == expected
    end

    @testset "NonlinearFunction" begin
        f_oop(u, p) = p[1] .* u
        f_iip(du, u, p) = (du .= p[1] .* u)

        nl_oop = NonlinearFunction{false}(f_oop)
        nl_iip = NonlinearFunction{true}(f_iip)

        expected = p[1] .* u

        # Standard calls
        @test nl_oop(u, p) == expected
        du = zeros(3)
        nl_iip(du, u, p)
        @test du == expected

        # Cross calls: IIP function with OOP style
        @test nl_iip(u, p) == expected

        # Cross calls: OOP function with IIP style
        du = zeros(3)
        nl_oop(du, u, p)
        @test du == expected
    end

    @testset "DynamicalODEFunction preserves standard behavior" begin
        # DynamicalODEFunction has internal ODEFunctions with non-standard signatures
        # (v, u, p, t) instead of (u, p, t). We need to ensure these still work correctly.
        dode_f1_op(v, u, p, t) = -u
        dode_f2_op(v, u, p, t) = p[1] .* v
        dode_f1_ip(dv, v, u, p, t) = (dv .= -u)
        dode_f2_ip(dv, v, u, p, t) = (dv .= p[1] .* v)

        dodefun_oop = DynamicalODEFunction{false}(dode_f1_op, dode_f2_op)
        dodefun_iip = DynamicalODEFunction{true}(dode_f1_ip, dode_f2_ip)

        v = [4.0, 5.0, 6.0]
        up = ArrayPartition(v, u)
        expected = ArrayPartition(dode_f1_op(v, u, p, t), dode_f2_op(v, u, p, t))

        # OOP call
        result = dodefun_oop(up, p, t)
        @test result.x[1] == expected.x[1]
        @test result.x[2] == expected.x[2]

        # IIP call
        dup = ArrayPartition(zeros(3), zeros(3))
        dodefun_iip(dup, up, p, t)
        @test dup.x[1] == expected.x[1]
        @test dup.x[2] == expected.x[2]
    end
end
