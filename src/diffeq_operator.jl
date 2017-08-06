### AbstractDiffEqOperator Interface

#=
1. Function call and multiplication: L(t,u,du) for inplace and du = L(t,u) for
   out-of-place, meaning L*u and A_mul_B!.
2. If the operator is not a constant, update it with (t,u). A mutating form, i.e.
   update_coefficients!(A,t,u) that changes the internal coefficients, and a
   out-of-place form B = update_coefficients(A,t,u).
3. is_constant(A) trait for whether the operator is constant or not.
=#

### AbstractDiffEqLinearOperator Interface

#=
1. AbstractDiffEqLinearOperator <: AbstractDiffEqOperator
2. Invariance under multiplication by a scalar
4. is_constant(A) trait for whether the operator is constant or not.
5. diagonal, symmetric, etc traits from LinearMaps.jl.
6. Optional: expm(A). Required for simple exponential integration.
7. Optional: expmv(A,t,u) = expm(t*A)*u and expmv!(v,A::DiffEqOperator,t,u)
   Required for sparse-saving exponential integration.
8. Optional: factorizations. A_ldiv_B, factorize et. al. This is only required
   for algorithms which use the factorization of the operator (Crank-Nicholson),
   and only for when the default linear solve is used.
=#

"""
AffineDiffEqOperator{T} <: AbstractDiffEqOperator{T}

`Ex: (A₁(t) + ... + Aₙ(t))*u + B₁(t) + ... + Bₙ(t)`

AffineDiffEqOperator{T}(As,Bs,u_cache=nothing)

Takes in two tuples for split Affine DiffEqs

1. update_coefficients! works by updating the coefficients of the component
   operators.
2. Function calls L(t,u) and L(t,u,du) are fallbacks interpretted in this form.
   This will allow them to work directly in the nonlinear ODE solvers without
   modification.
3. f(t,u,du) is only allowed if a u_cache is given
4. B(t) can be Union{Number,AbstractArray}, in which case they are constants.
   Otherwise they are interpreted they are functions v=B(t) and B(v,t)

Solvers will see this operator from integrator.f and can interpret it by
checking the internals of As and Bs. For example, it can check is_constant(As[1])
etc.
"""
struct AffineDiffEqOperator{T,T1,T2,U} <: AbstractDiffEqOperator{T}
    As::T1
    Bs::T2
    u_cache::U
    function AffineDiffEqOperator{T}(As,Bs,u_cache=nothing) where T
        new{T,typeof(As),typeof(Bs),typeof(u_cache)}(As,Bs,u_cache)
    end
end


function (L::AffineDiffEqOperator)(t,u)
    tmp = sum((update_coefficients(A,t,u); A*u for A in L.As))
    tmp2 = sum((typeof(B) <: Union{Number,AbstractArray} ? B : B(t) for B in L.Bs))
    tmp + tmp2
end

function (L::AffineDiffEqOperator)(t,u,du)
    update_coefficients!(L,t,u)
    L.u_cache == nothing && error("Can only use inplace AffineDiffEqOperator if u_cache is given.")
    u_cache = L.u_cache
    fill!(du,zero(du))
    # TODO: Make type-stable via recursion
    for A in L.As
        A_mul_B!(u_cache,A,u)
        du .+= u_cache
    end
    for B in L.Bs
        if typeof(B) <: Union{Number,AbstractArray}
            du .+= B
        else
            B(u_cache,t)
            du .+= u_cache
        end
    end
end

function update_coefficients!(L::AffineDiffEqOperator,t,u)
    # TODO: Make type-stable via recursion
    for A in L.As; update_coefficients!(A,t,u); end
    for B in L.Bs; update_coefficients!(B,t,u); end
end

"""
DiffEqScalar Interface

DiffEqScalar(func,coeff=1.0)

This is a function with a coefficient.

α(t) returns a new DiffEqScalar with an updated coefficient.
"""
struct DiffEqScalar{F,T}
    func::F
    coeff::T
    DiffEqScalar{T}(func) where T = new{typeof(func),T}(func,one(T))
    DiffEqScalar{F,T}(func,coeff) where {F,T} = new{F,T}(func,coeff)
end
DiffEqScalar(func,coeff=1.0) = DiffEqScalar{typeof(func),typeof(coeff)}(func,coeff)
function (α::DiffEqScalar)(t)
    DiffEqScalar(α.func,α.func(t))
end
Base.:*(α::Number,B::DiffEqScalar) = DiffEqScalar(B.func,B.coeff*α)
Base.:*(B::DiffEqScalar,α::Number) = DiffEqScalar(B.func,B.coeff*α)

# Inherits the standard assumptions of an AbstractLinearMap
# Extra standard assumptions
is_constant(L::AbstractDiffEqLinearOperator) = true
update_coefficients!(L,t,u) = nothing
update_coefficients(L,t,u) = L

# Generic fallbacks
Base.expm(L::AbstractDiffEqLinearOperator,t) = expm(t*L)
expmv(L::AbstractDiffEqLinearOperator,t,u) = expm(t,L)*u
expmv!(v,L::AbstractDiffEqLinearOperator,t,u) = A_mul_B!(v,expm(t,L),u)

### AbstractDiffEqLinearOperator defined by an array and update functions
struct DiffEqArrayOperator{T,Arr<:AbstractMatrix{T},Sca,F} <: AbstractDiffEqLinearOperator{T}
    A::Arr
    α::Sca
    _isreal::Bool
    _issymmetric::Bool
    _ishermitian::Bool
    _isposdef::Bool
    update_func::F
end
DEFAULT_UPDATE_FUNC = (L,t,u)->nothing
function DiffEqArrayOperator(A::AbstractMatrix{T},α=1.0,
                             update_func = DEFAULT_UPDATE_FUNC) where T
    if !(typeof(α) <: Number)
        _α = DiffEqScalar(nothing,α)
    elseif !(typeof(α) <: DiffEqScalar) # Assume it's some kind of function
        _α = DiffEqScalar(α,1.0)
    else # Must be a DiffEqScalar already
        _α = α
    end
    DiffEqArrayOperator{T,typeof(A),typeof(_α),
    typeof(update_func)}(
    A,_α,isreal(A),issymmetric(A),ishermitian(A),
    isposdef(A),update_func)
end

Base.isreal(L::DiffEqArrayOperator) = L._isreal
Base.issymmetric(L::DiffEqArrayOperator) = L._issymmetric
Base.ishermitian(L::DiffEqArrayOperator) = L._ishermitian
Base.isposdef(L::DiffEqArrayOperator) = L._isposdef
is_constant(L::DiffEqArrayOperator) = L.update_func == DEFAULT_UPDATE_FUNC

update_coefficients!(L::DiffEqArrayOperator,t,u) = (L.update_func(L.A,t,u); L.α = L.α(t); nothing)
update_coefficients(L::DiffEqArrayOperator,t,u)  = (L.update_func(L.A,t,u); L.α = L.α(t); L)

function (L::DiffEqArrayOperator)(t,u)
  update_coefficients!(L,t,u)
  L*u
end
function (L::DiffEqArrayOperator)(t,u,du)
  update_coefficients!(L,t,u)
  A_mul_B!(du,L,u)
end

### Forward some extra operations
function Base.:*(α::Number,L::DiffEqArrayOperator)
    DiffEqArrayOperator(L.A,DiffEqScalar(L.func,L.coeff*α),L.update_func)
end
Base.:*(L::DiffEqArrayOperator,α::Number) = α*L
Base.:*(L::DiffEqArrayOperator,b::AbstractArray) = L.α.coeff*L.A*b
function Base.A_mul_B!(v,L::DiffEqArrayOperator,b::AbstractArray)
    A_mul_B!(v,L.A,b)
    scale!(b,L.α.coeff)
end
Base.expm(L::DiffEqArrayOperator) = expm(L.α.coeff*L.A)
Base.size(L::DiffEqArrayOperator) = size(L.A)

function Base.A_ldiv_B!(x,L::DiffEqArrayOperator, b::AbstractArray)
    A_ldiv_B!(x,L.A,b)
    scale!(x,inv(L.α.coeff))
end

"""
FactorizedDiffEqArrayOperator{T,I}

A helper function for holding factorized version of the DiffEqArrayOperator
"""
struct FactorizedDiffEqArrayOperator{T,I}
    A::T
    inv_coeff::I
end

Base.factorize(L::DiffEqArrayOperator)         = FactorizedDiffEqArrayOperator(factorize(L.A),inv(L.α.coeff))
Base.lufact(L::DiffEqArrayOperator,args...)    = FactorizedDiffEqArrayOperator(lufact(L.A,args...),inv(L.α.coeff))
Base.lufact!(L::DiffEqArrayOperator,args...)   = FactorizedDiffEqArrayOperator(lufact!(L.A,args...),inv(L.α.coeff))
Base.qrfact(L::DiffEqArrayOperator,args...)    = FactorizedDiffEqArrayOperator(qrfact(L.A,args...),inv(L.α.coeff))
Base.qrfact!(L::DiffEqArrayOperator,args...)   = FactorizedDiffEqArrayOperator(qrfact!(L.A,args...),inv(L.α.coeff))
Base.cholfact(L::DiffEqArrayOperator,args...)  = FactorizedDiffEqArrayOperator(cholfact(L.A,args...),inv(L.α.coeff))
Base.cholfact!(L::DiffEqArrayOperator,args...) = FactorizedDiffEqArrayOperator(cholfact!(L.A,args...),inv(L.α.coeff))
Base.ldltfact(L::DiffEqArrayOperator,args...)  = FactorizedDiffEqArrayOperator(ldltfact(L.A,args...),inv(L.α.coeff))
Base.ldltfact!(L::DiffEqArrayOperator,args...) = FactorizedDiffEqArrayOperator(ldltfact!(L.A,args...),inv(L.α.coeff))
Base.bkfact(L::DiffEqArrayOperator,args...)    = FactorizedDiffEqArrayOperator(bkfact(L.A,args...),inv(L.α.coeff))
Base.bkfact!(L::DiffEqArrayOperator,args...)   = FactorizedDiffEqArrayOperator(bkfact!(L.A,args...),inv(L.α.coeff))
Base.lqfact(L::DiffEqArrayOperator,args...)    = FactorizedDiffEqArrayOperator(lqfact(L.A,args...),inv(L.α.coeff))
Base.lqfact!(L::DiffEqArrayOperator,args...)   = FactorizedDiffEqArrayOperator(lqfact!(L.A,args...),inv(L.α.coeff))
Base.svdfact(L::DiffEqArrayOperator,args...)   = FactorizedDiffEqArrayOperator(svdfact(L.A,args...),inv(L.α.coeff))
Base.svdfact!(L::DiffEqArrayOperator,args...)  = FactorizedDiffEqArrayOperator(svdfact!(L.A,args...),inv(L.α.coeff))

function Base.A_ldiv_B!(x,L::FactorizedDiffEqArrayOperator, b::AbstractArray)
    A_ldiv_B!(x,L.A,b)
    scale!(x,inv(L.inv_coeff))
end
function Base.:\(L::FactorizedDiffEqArrayOperator, b::AbstractArray)
    (L.A \ b) * L.inv_coeff
end
