using DiffEqBase, ForwardDiff, Test

u0 = 2.0
p = 2.0
t0 = 1.0

@test DiffEqBase.promote_u0(u0,p,t0) isa Float64
@test DiffEqBase.promote_u0(u0,p,t0) == 2.0

p_possibilities = [ForwardDiff.Dual(2.0),(ForwardDiff.Dual(2.0),2.0),
    [ForwardDiff.Dual(2.0)],([ForwardDiff.Dual(2.0)],2.0),
    (2.0,ForwardDiff.Dual(2.0)),(;x=2.0,y=ForwardDiff.Dual(2.0)),
    (;x=2.0,y=[ForwardDiff.Dual(2.0)]),(;x=2.0,y=[[ForwardDiff.Dual(2.0)]])
]

for p in p_possibilities
    @test DiffEqBase.anyeltypedual(p) <: ForwardDiff.Dual
    u0 = 2.0
    @test DiffEqBase.promote_u0(u0,p,t0) isa ForwardDiff.Dual
    u0 = ForwardDiff.Dual(2.0)
    @test DiffEqBase.promote_u0(u0,p,t0) isa ForwardDiff.Dual
    @inferred DiffEqBase.anyeltypedual(p)
end

p_possibilities_uninferrred = [
    Dict(:x => 2.0, :y => ForwardDiff.Dual(2.0)),
    Dict(:x => 2.0, :y => [ForwardDiff.Dual(2.0)]),
    Dict(:x => 2.0, :y => [(;x=(ForwardDiff.Dual(2.0),2.0),y=2.0)]),
]

for p in p_possibilities_uninferrred
    @test DiffEqBase.anyeltypedual(p) <: ForwardDiff.Dual
    u0 = 2.0
    @test DiffEqBase.promote_u0(u0,p,t0) isa ForwardDiff.Dual
    u0 = ForwardDiff.Dual(2.0)
    @test DiffEqBase.promote_u0(u0,p,t0) isa ForwardDiff.Dual
end