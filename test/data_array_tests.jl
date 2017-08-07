using DiffEqBase, Base.Test

type VectorType{T} <: DEDataArray{T}
    x::Vector{T}
    f::T
end

type MatrixType{T,S} <: DEDataArray{T}
    x::Matrix{T}
    f::S
end

a = VectorType{Float64}([0.0; 1.0], 2.0)
b = MatrixType{Int,Float64}([1  2; 4  3], 3.0)

# basic methods of AbstractArray interface
@test eltype(a) == Float64 && eltype(b) == Int
@test length(a) == 2 && size(a) == (2,) && eachindex(a) == Base.OneTo(2)
@test length(b) == 4 && size(b) == (4,) && eachindex(b) == Base.OneTo(4)
@test first(a) == 0.0 && last(a) == 1.0
@test first(b) == 1 && last(b) == 3

# simple broadcasts
@test a .+ [1  2; 4  3] == [1.0  2.0; 5.0  4.0]
@test [0.0; 1.0] .+ b == [1.0  2.0; 5.0  4.0]
@test [1//2] .* a == [0.0; 0.5]
@test b ./ 2 ==  [0.5  1.0; 2.0  1.5]
@test_throws MethodError a .+ b

a2 = similar(a)
@test a.x == a2.x && a.f == a2.f

# broadcast assignments
a2.f = 3
a .= a .+ a2
@test a == VectorType{Float64}([0.0; 2.0], 2.0)

b .= b.^2
@test b == MatrixType{Int,Float64}([1  4; 16  9], 3.0)
