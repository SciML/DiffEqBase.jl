# Note: This must be the first file executed in the tests.
#
# The following structure is used to test the problem reported at issue #507. We
# need to define it here because, since `recursivecopy!` and `copy_fields!` are
# generated functions, then the `ArrayInterface.isimmutable` method must be
# defined before the first `using DiffEqBase` in a Julia session.

struct Quaternion{T} <: AbstractVector{T}
    q0::T
    q1::T
    q2::T
    q3::T
end

using ArrayInterface
ArrayInterface.ismutable(::Type{<:Quaternion}) = false
Base.size(::Quaternion) = 4

using DiffEqBase, RecursiveArrayTools, Test

mutable struct VectorType{T} <: DEDataVector{T}
    x::Vector{T}
    f::T
end

mutable struct MatrixType{T,S} <: DEDataMatrix{T}
    f::S
    x::Matrix{T}
end

A = [0.0; 1.0]
B = [1.0  2; 4  3]

a = VectorType{Float64}(copy(A), 1.0)
b = MatrixType{Float64,Float64}(2.0, copy(B))

# basic methods of AbstractArray interface
@test eltype(a) == Float64 && eltype(b) == Float64

# size
@test length(a) == 2 && length(b) == 4
@test size(a) == (2,) && size(b) == (2,2)

# iteration
@test first(a) == 0.0 && first(b) == 1
@test last(a) == 1.0 && last(b) == 3
for (i, (ai, Ai)) in enumerate(zip(a, A))
    @test ai == Ai
end
for (i, (bi, Bi)) in enumerate(zip(b, B))
    @test bi == Bi
end

# indexing
@test eachindex(a) == Base.LinearIndices(A) && eachindex(b) == vec(Base.LinearIndices(B))
@test a[2] == a[end] == 1.0
@test b[3] == b[1,2] == 2 && b[:, 1] == [1; 4]

a[1] = 3; b[2,1] = 1
@test a[1] == 3.0 && b[2] == 1

a[:] = A; b[:, 1] = B[1:2]
@test a.x == A && b.x == B

# simple broadcasts
@test [1//2] .* a == [0.0; 0.5]
@test b ./ 2 ==  [0.5  1.0; 2.0  1.5]
@test b .+ a  == [1.0  2.0; 5.0  4.0]
@test_broken a .+ b  == [1.0  2.0; 5.0  4.0] # Doesn't find the largest

# similar data arrays
a2 = similar(a); b2 = similar(b, (1,4)); b3 = similar(b, Float64, (1, 4))
@test typeof(a2.x) == Vector{Float64} && a2.f == 1.0
@test typeof(b2.x) == Matrix{Float64} && b2.f == 2.0
@test typeof(b3.x) == Matrix{Float64} && b3.f == 2.0

# copy all fields of data arrays
recursivecopy!(a2, a)
recursivecopy!(b2, b)
@test a2.x == A && vec(b2.x) == vec(B)
recursivecopy!(b3, b)

# copy all fields except of field x
a2.f = 3.0; b2.f = -1; b3.f == 0
DiffEqBase.copy_fields!(a, a2)
DiffEqBase.copy_fields!(b, b2)
@test a.f == 3.0 && b.f == -1.0
DiffEqBase.copy_fields!(b, b3)

# create data array with field x replaced by new array
a3 = DiffEqBase.copy_fields([1.0; 0.0], a)
@test a3 == VectorType{Float64}([1; 0], 3.0)

# broadcast assignments
a.f = 0.0
a .= a2 .+ a3
@test a == VectorType{Float64}([1.0; 1.0], 0.0)
@test a == a2 .+ a3
@test (a2 .+ a3) isa VectorType

old_b = copy(b)
b .= b .^ 2 .+ a3
@test b == MatrixType{Int,Float64}(-1.0, [2  5; 16  9])
@test b == old_b .^ 2 .+ a3
@test (b .^ 2 .+ a3) isa MatrixType

using StaticArrays

# Test ability to use MVectors in the solvers
mutable struct SimWorkspace{T} <: DEDataVector{T}
  x::MVector{2,T}
  a::T
end
s0   = SimWorkspace{Float64}(MVector{2,Float64}(1.0,4.0),1.)
similar(s0,Float64,size(s0))
s1   = SimWorkspace{Float64}(MVector{2,Float64}(2.0,1.0),1.)
s0 .+ s1 == SimWorkspace{Float64}(MVector{2,Float64}(3.0,5.0),1.)

mutable struct SimWorkspace2{T} <: DEDataVector{T}
  x::SVector{2,T}
  a::T
end
s0   = SimWorkspace2{Float64}(SVector{2,Float64}(1.0,4.0),1.)
s1   = SimWorkspace2{Float64}(SVector{2,Float64}(2.0,1.0),1.)
s0 .+ s1 == SimWorkspace2{Float64}(SVector{2,Float64}(3.0,5.0),1.)

# Test `recursivecopy!` in immutable structures derived from `AbstractArrays`.
# See issue #507.
mutable struct SimWorkspace3{T} <: DEDataVector{T}
    x::Vector{T}
    q::Quaternion{T}
end

a = SimWorkspace3([1.0,2.0,3.0], Quaternion(cosd(15), 0.0, 0.0, sind(15)))
b = SimWorkspace3([0.0,0.0,0.0], Quaternion(1.0, 0.0, 0.0, 0.0))

recursivecopy!(b,a)
@test b.x == a.x
@test b.q.q0 == a.q.q0
@test b.q.q1 == a.q.q1
@test b.q.q2 == a.q.q2
@test b.q.q3 == a.q.q3

a = SimWorkspace3([1.0,2.0,3.0], Quaternion(cosd(15), 0.0, 0.0, sind(15)))
b = SimWorkspace3([0.0,0.0,0.0], Quaternion(1.0, 0.0, 0.0, 0.0))

DiffEqBase.copy_fields!(b, a)
@test b.x == [0.0,0.0,0.0]
@test b.q.q0 == a.q.q0
@test b.q.q1 == a.q.q1
@test b.q.q2 == a.q.q2
@test b.q.q3 == a.q.q3

