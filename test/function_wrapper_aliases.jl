using DiffEqBase, ForwardDiff, Test
using DiffEqBase: Void, FunctionWrappersWrappers, OrdinaryDiffEqTag,
    wrapfun_iip, wrapfun_iip_simple, AnyFunctionWrapper
using FunctionWrappers: FunctionWrapper

# Get the ForwardDiff extension module
const FDExt = Base.get_extension(DiffEqBase, :DiffEqBaseForwardDiffExt)

@testset "DEIIPFunctionWrapper struct (no ForwardDiff)" begin
    ff = (du, u, p, t) -> (du .= u; nothing)
    du = zeros(3)
    u = ones(3)
    p = [1.0, 2.0]
    t = 0.0

    # wrapfun_iip_simple should produce a DEIIPFunctionWrapper
    wrapped = wrapfun_iip_simple(ff, du, u, p, t)
    @test wrapped isa DiffEqBase.DEIIPFunctionWrapper
    @test wrapped isa DiffEqBase.DEIIPFunctionWrapper{
        Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64}

    # VF64 alias should match
    @test wrapped isa DiffEqBase.DEIIPFunctionWrapperVF64{Vector{Float64}}

    # Should match AnyFunctionWrapper union
    @test wrapped isa AnyFunctionWrapper

    # The wrapped function should be callable and produce correct results
    du_test = zeros(3)
    wrapped(du_test, u, p, t)
    @test du_test == u

    # isfunctionwrapper should return true
    @test SciMLBase.isfunctionwrapper(wrapped)

    # Stack trace type string should be short
    type_str = string(typeof(wrapped))
    @test occursin("DEIIPFunctionWrapper", type_str)
    @test !occursin("FunctionWrappersWrapper", type_str)
end

@testset "DEIIPFunctionWrapperVF64 with NullParameters" begin
    ff = (du, u, p, t) -> (du .= u; nothing)
    du = zeros(3)
    u = ones(3)
    p = SciMLBase.NullParameters()
    t = 0.0

    wrapped = wrapfun_iip_simple(ff, du, u, p, t)
    @test wrapped isa DiffEqBase.DEIIPFunctionWrapper
    @test wrapped isa DiffEqBase.DEIIPFunctionWrapperVF64{SciMLBase.NullParameters}
end

@testset "ODEDualTag and ODEDualType (ForwardDiff extension)" begin
    @test FDExt.ODEDualTag === ForwardDiff.Tag{OrdinaryDiffEqTag, Float64}
    @test FDExt.ODEDualType === ForwardDiff.Dual{
        ForwardDiff.Tag{OrdinaryDiffEqTag, Float64}, Float64, 1}
end

@testset "DEIIPFunctionWrapperForwardDiff struct" begin
    ff = (du, u, p, t) -> (du .= u; nothing)
    du = zeros(3)
    u = ones(3)
    p = [1.0, 2.0]
    t = 0.0

    # wrapfun_iip with ForwardDiff loaded should produce a DEIIPFunctionWrapperForwardDiff
    wrapped = wrapfun_iip(ff, (du, u, p, t))
    @test wrapped isa DiffEqBase.DEIIPFunctionWrapperForwardDiff

    # VF64 alias should match
    @test wrapped isa FDExt.DEIIPFunctionWrapperForwardDiffVF64{Vector{Float64}}

    # Should match AnyFunctionWrapper union
    @test wrapped isa AnyFunctionWrapper

    # The wrapped function should be callable
    du_test = zeros(3)
    wrapped(du_test, u, p, t)
    @test du_test == u

    # isfunctionwrapper should return true
    @test SciMLBase.isfunctionwrapper(wrapped)

    # Stack trace type string should NOT contain FunctionWrappersWrapper
    type_str = string(typeof(wrapped))
    @test occursin("DEIIPFunctionWrapperForwardDiff", type_str)
    @test !occursin("FunctionWrappersWrapper", type_str)
end

@testset "DEIIPFunctionWrapperForwardDiffVF64 with NullParameters" begin
    ff = (du, u, p, t) -> (du .= u; nothing)

    # Default wrapfun_iip (no args) produces the 7-wrapper variant, not 4-wrapper
    wrapped_default = wrapfun_iip(ff)
    # The default 7-wrapper has a different structure (7 entries, not 4),
    # so it should NOT be a DEIIPFunctionWrapperForwardDiff
    @test !(wrapped_default isa DiffEqBase.DEIIPFunctionWrapperForwardDiff)
    # But it IS still an AnyFunctionWrapper (raw FunctionWrappersWrapper)
    @test wrapped_default isa AnyFunctionWrapper

    # With explicit 4-tuple args and NullParameters, it should match
    du = zeros(3)
    u = ones(3)
    p = SciMLBase.NullParameters()
    t = 0.0
    wrapped = wrapfun_iip(ff, (du, u, p, t))
    @test wrapped isa DiffEqBase.DEIIPFunctionWrapperForwardDiff
    @test wrapped isa FDExt.DEIIPFunctionWrapperForwardDiffVF64{SciMLBase.NullParameters}
end

@testset "wrapfun_iip_simple does not change behavior with ForwardDiff loaded" begin
    ff = (du, u, p, t) -> (du .= u; nothing)
    du = zeros(3)
    u = ones(3)
    p = [1.0, 2.0]
    t = 0.0

    # wrapfun_iip_simple should ALWAYS produce a DEIIPFunctionWrapper, even with ForwardDiff loaded
    wrapped_simple = wrapfun_iip_simple(ff, du, u, p, t)
    @test wrapped_simple isa DiffEqBase.DEIIPFunctionWrapper

    # It should NOT match the ForwardDiff 4-wrapper struct
    @test !(wrapped_simple isa DiffEqBase.DEIIPFunctionWrapperForwardDiff)

    # wrapfun_iip should produce a DEIIPFunctionWrapperForwardDiff (ForwardDiff ext wraps it)
    wrapped_fd = wrapfun_iip(ff, (du, u, p, t))
    @test wrapped_fd isa DiffEqBase.DEIIPFunctionWrapperForwardDiff
    @test !(wrapped_fd isa DiffEqBase.DEIIPFunctionWrapper)
end
