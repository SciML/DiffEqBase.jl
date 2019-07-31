using DiffEqBase, Test

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


function test_inplace(du, expected, f::Function, args...)
    """Test the in-place version of a function."""
    fill!(du, NaN)
    f(du, args...)
    @test du == expected
end

# Allocate du automatically based on type of expected result
test_inplace(expected, f::Function, args...) = test_inplace(similar(expected), expected, f, args...)

function test_iop(expected, f_op::Function, f_ip::Function, args...)
    """Test in- and out-of-place version of function both match expected value."""
    @test f_op(args...) == expected
    test_inplace(expected, f_ip, args...)
end


@iop_def f(u,p,t) = p[1] .* u
u = [1.0, 2.0, 3.0]
p = [2.0]
t = 0.0


# ODEFunction
odefun = ODEFunction{false}(f_op)
odefun_ip = ODEFunction{true}(f_ip)
expected = f_op(u, p, t)
test_iop(expected, odefun, odefun_ip, u, p, t)

# DiscreteFunction
dfun = DiscreteFunction{false}(f_op)
dfun_ip = DiscreteFunction{true}(f_ip)
test_iop(f_op(u, p, t), dfun, dfun_ip, u, p, t)


# Type stability
f_analytic(u,p,t) = u
jac = (u,p,t) -> 1
@inferred ODEFunction{false}(f_op,jac=jac)
@inferred DiscreteFunction{false}(f_op,analytic=f_analytic)
