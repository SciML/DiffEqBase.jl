using DiffEqBase, Test

mutable struct VectorType{T} <: DEDataVector{T}
    x::Vector{T}
    f::T
end

mutable struct MatrixType{T,S} <: DEDataMatrix{T}
    f::S
    x::Matrix{T}
end

A = [0.0; 1.0]
B = [1  2; 4  3]

a = VectorType{Float64}(copy(A), 1.0)
b = MatrixType{Int,Float64}(2.0, copy(B))

# basic methods of AbstractArray interface
@test eltype(a) == Float64 && eltype(b) == Int

# size
@test length(a) == 2 && length(b) == 4
@test size(a) == (2,) && size(b) == (2,2)

# iteration
@test first(a) == 0.0 && first(b) == 1
@test last(a) == 1.0 && last(b) == 3

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
@test a .+ b == [1.0  2.0; 5.0  4.0]

# similar data arrays
a2 = similar(a); b2 = similar(b, (1,4)); b3 = similar(b, Float64, (1, 4))
@test typeof(a2.x) == Vector{Float64} && a2.f == 1.0
@test typeof(b2.x) == Matrix{Int} && b2.f == 2.0
@test typeof(b3.x) == Matrix{Float64} && b3.f == 2.0

# copy all fields of data arrays
DiffEqBase.recursivecopy!(a2, a)
DiffEqBase.recursivecopy!(b2, b)
@test a2.x == A && vec(b2.x) == vec(B)
DiffEqBase.recursivecopy!(b3, b)

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

b .= b .^ 2 .+ a3
@test b == MatrixType{Int,Float64}(-1.0, [2  5; 16  9])

using StaticArrays

# Test ability to use MVectors in the solvers
mutable struct SimWorkspace{T} <: DEDataVector{T}
  x::MVector{2,T}
  a::T
end
s0   = SimWorkspace{Float64}(MVector{2,Float64}(0,0),1.)
similar(s0,Float64,size(s0))
