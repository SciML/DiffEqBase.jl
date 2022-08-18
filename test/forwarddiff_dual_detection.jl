using DiffEqBase, ForwardDiff, Test, InteractiveUtils
using Plots

u0 = 2.0
p = 2.0
t0 = 1.0

@test DiffEqBase.promote_u0(u0, p, t0) isa Float64
@test DiffEqBase.promote_u0(u0, p, t0) == 2.0

struct MyStruct{T, T2} <: Number
    x::T
    y::T2
end

struct MyStruct2{T, T2}
    x::T
    y::T2
    MyStruct2(x) = new{typeof(x), Any}(x)
end

struct MyStruct3{T, T2}
    x::T
    y::T2
    MyStruct3(x) = new{typeof(x), Float64}(x)
end

module Mod end

p_possibilities = [ForwardDiff.Dual(2.0), (ForwardDiff.Dual(2.0), 2.0),
    [ForwardDiff.Dual(2.0)], ([ForwardDiff.Dual(2.0)], 2.0),
    (2.0, ForwardDiff.Dual(2.0)), (; x = 2.0, y = ForwardDiff.Dual(2.0)),
    (; x = 2.0, y = [ForwardDiff.Dual(2.0)]), (; x = 2.0, y = [[ForwardDiff.Dual(2.0)]]),
    Set([2.0, ForwardDiff.Dual(2.0)]), (SciMLBase.NullParameters(), ForwardDiff.Dual(2.0)),
    ((), ForwardDiff.Dual(2.0)), ForwardDiff.Dual{Nothing}(ForwardDiff.Dual{MyStruct}(2.0)),
    (plot(), ForwardDiff.Dual(2.0)),
]

for p in p_possibilities
    @show p
    @test DiffEqBase.anyeltypedual(p) <: ForwardDiff.Dual
    u0 = 2.0
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    u0 = ForwardDiff.Dual(2.0)
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    @inferred DiffEqBase.anyeltypedual(p)
end

higher_order_p_possibilities = [ForwardDiff.Dual{Nothing}(ForwardDiff.Dual{MyStruct}(2.0)),
    (ForwardDiff.Dual{Nothing}(ForwardDiff.Dual{MyStruct}(2.0)),
     SciMLBase.NullParameters()),
    (ForwardDiff.Dual{Nothing}(ForwardDiff.Dual{MyStruct}(2.0)),
     ForwardDiff.Dual{Nothing}(2.0)),
    (ForwardDiff.Dual{Nothing}(2.0),
     ForwardDiff.Dual{Nothing}(ForwardDiff.Dual{MyStruct}(2.0))),
]

for p in higher_order_p_possibilities
    @show p
    @test DiffEqBase.anyeltypedual(p) <: ForwardDiff.Dual
    @test DiffEqBase.anyeltypedual(p) <:
          ForwardDiff.Dual{Nothing, ForwardDiff.Dual{MyStruct, Float64, 0}, 0}
    u0 = 2.0
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    u0 = ForwardDiff.Dual(2.0)
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    @inferred DiffEqBase.anyeltypedual(p)
end

p_possibilities17 = [
    MyStruct(2.0, ForwardDiff.Dual(2.0)), [MyStruct(2.0, ForwardDiff.Dual(2.0))],
    [MyStruct(2.0, [2.0, ForwardDiff.Dual(2.0)])],
    [MyStruct(2.0, (2.0, ForwardDiff.Dual(2.0)))],
    ((;), ForwardDiff.Dual(2.0)), MyStruct3(ForwardDiff.Dual(2.0)),
    (Mod, ForwardDiff.Dual(2.0)), (() -> 2.0, ForwardDiff.Dual(2.0)),
    (Base.pointer([2.0]), ForwardDiff.Dual(2.0)),
]
VERSION >= v"1.7" &&
    push!(p_possibilities17, Returns((a = 2, b = 1.3, c = ForwardDiff.Dual(2.0f0))))

for p in p_possibilities17
    @show p
    @test DiffEqBase.anyeltypedual(p) <: ForwardDiff.Dual
    u0 = 2.0
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    u0 = ForwardDiff.Dual(2.0)
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual

    if VERSION >= v"1.7"
        # v1.6 does not infer `getproperty` mapping
        @inferred DiffEqBase.anyeltypedual(p)
        ci = InteractiveUtils.@code_typed DiffEqBase.anyeltypedual(p)
        @show filter(!=(Expr(:code_coverage_effect)), ci.first.code)
        @test count(x -> (x != (Expr(:code_coverage_effect))) &&
                        (x != GlobalRef(DiffEqBase, :Any)), ci.first.code) == 1
    end
end

p_possibilities_uninferrred = [
    Dict(:x => 2.0, :y => ForwardDiff.Dual(2.0)),
    Dict(:x => 2.0, :y => [ForwardDiff.Dual(2.0)]),
    Dict(:x => 2.0, :y => [(; x = (ForwardDiff.Dual(2.0), 2.0), y = 2.0)]),
    Dict(:x => 2.0, :y => [(; x = [MyStruct(2.0, [2.0, ForwardDiff.Dual(2.0)])], y = 2.0)]),
    [MyStruct("2", [2.0, ForwardDiff.Dual(2.0)])],
    Dict(:x => [MyStruct("2", [2.0, MyStruct(ForwardDiff.Dual(2.0), 2.0)])],
         :y => ForwardDiff.Dual{MyStruct}(2.0)),
    ((Dict(:x => nothing)), ForwardDiff.Dual(2.0)),
    MyStruct2(ForwardDiff.Dual(2.0)),
    [MyStruct2(ForwardDiff.Dual(2.0)), 2.0],

    # Vectors of non-number types won't infer
    [MyStruct(2.0, ForwardDiff.Dual(2.0))],
    (; x = 2.0, y = [[MyStruct3(ForwardDiff.Dual(2.0))]]),
    (; x = Vector{Float64}(undef, 2), y = [[MyStruct3(ForwardDiff.Dual(2.0))]]),
    (; x = Matrix{Any}(undef, 2, 2), y = [[MyStruct3(ForwardDiff.Dual(2.0))]]),
]

for p in p_possibilities_uninferrred
    @show p
    @test DiffEqBase.anyeltypedual(p) <: ForwardDiff.Dual
    u0 = 2.0
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    u0 = ForwardDiff.Dual(2.0)
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
end

p_possibilities_missed = [
    Set([2.0, "s", ForwardDiff.Dual(2.0)]),
    Set([2.0, ForwardDiff.Dual(2.0), SciMLBase.NullParameters()]),
    Set([Matrix{Float64}(undef, 2, 2), ForwardDiff.Dual(2.0)]),
]

for p in p_possibilities_missed
    @show p
    @test_broken DiffEqBase.anyeltypedual(p) <: ForwardDiff.Dual
    u0 = 2.0
    @test_broken DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    u0 = ForwardDiff.Dual(2.0)
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
end

p_possibilities_notdual = [
    (), (;), [2.0], [2.0, 2], [2.0, (2.0)], [2.0, MyStruct(2.0, 2.0f0)],
]

for p in p_possibilities_notdual
    @show p
    @test !(DiffEqBase.anyeltypedual(p) <: ForwardDiff.Dual)
    u0 = 2.0
    @test !(DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual)
    u0 = ForwardDiff.Dual(2.0)
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    @inferred DiffEqBase.anyeltypedual(p)
end

p_possibilities_notdual_uninferred = [
    [],

    # Undefs cause inference loss
    [2.0, MyStruct3(2.0)], [2.0, MyStruct2(2.0)], [2.0, MyStruct2(2.0), []],
    [Dict(:x => 2, "y" => 5), MyStruct2(2.0)],

    # Dictionaries can have inference issues
    Dict(:x => 2, :y => 5), Dict(:x => 2, "y" => 5),
]

# Also check circular references
# https://github.com/SciML/DiffEqBase.jl/issues/784

x = Any[[1.0, 2.0]]
push!(x, x)
push!(p_possibilities_notdual_uninferred, x)

struct X
    x::Any
end
x = Any[[1.0, 2.0]]
push!(x, X(x))
push!(p_possibilities_notdual_uninferred, x)

mutable struct Y
    x::Any
end
x = Y(1)
x.x = x
push!(p_possibilities_notdual_uninferred, x)

for p in p_possibilities_notdual_uninferred
    @test !(DiffEqBase.anyeltypedual(p) <: ForwardDiff.Dual)
    u0 = 2.0
    @test !(DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual)
    u0 = ForwardDiff.Dual(2.0)
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
end

f(du, u, p, t) = du .= u
config = ForwardDiff.JacobianConfig(f, ones(5))

p_possibilities_configs = [
    (config, config), (config, 2.0), config, (; x = config, y = 2.0),
]

for p in p_possibilities_configs
    @show p
    @test !(DiffEqBase.anyeltypedual(p) <: ForwardDiff.Dual)
    u0 = 2.0
    @test !(DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual)
    u0 = ForwardDiff.Dual(2.0)
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
    @inferred DiffEqBase.anyeltypedual(p)
end

p_possibilities_configs_not_inferred = [
    [2.0, (2.0,), config], [2.0, config, MyStruct(2.0, 2.0f0)],
]

for p in p_possibilities_configs_not_inferred
    @show p
    @test !(DiffEqBase.anyeltypedual(p) <: ForwardDiff.Dual)
    u0 = 2.0
    @test !(DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual)
    u0 = ForwardDiff.Dual(2.0)
    @test DiffEqBase.promote_u0(u0, p, t0) isa ForwardDiff.Dual
end
