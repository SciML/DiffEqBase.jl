using LinearAlgebra
### AbstractDiffEqOperator Interface

#=
1. Function call and multiplication: L(u,p,t,du) for inplace and du = L(u,p,t) for
   out-of-place, meaning L*u and A_mul_B!.
2. If the operator is not a constant, update it with (u,p,t). A mutating form, i.e.
   update_coefficients!(A,u,p,t) that changes the internal coefficients, and a
   out-of-place form B = update_coefficients(A,u,p,t).
3. is_constant(A) trait for whether the operator is constant or not.
=#

Base.eltype(L::AbstractDiffEqOperator{T}) where T = T
update_coefficients!(L,u,p,t) = nothing
update_coefficients(L,u,p,t) = L

# Traits
is_constant(L::AbstractDiffEqOperator) = false
has_expmv!(L::AbstractDiffEqOperator) = false # expmv!(v, L, t, u)
has_expmv(L::AbstractDiffEqOperator) = false # v = exp(L, t, u)
has_exp(L::AbstractDiffEqOperator) = false # v = exp(L, t)*u
has_mul(L::AbstractDiffEqOperator) = true # du = L*u
has_mul!(L::AbstractDiffEqOperator) = false # A_mul_B!(du, L, u)
has_ldiv(L::AbstractDiffEqOperator) = false # du = L\u
has_ldiv!(L::AbstractDiffEqOperator) = false # ldiv!(du, L, u)

### AbstractDiffEqLinearOperator Interface

#=
1. AbstractDiffEqLinearOperator <: AbstractDiffEqOperator
2. Can absorb under multiplication by a scalar. In all algorithms things like
   dt*L show up all the time, so the linear operator must be able to absorb
   such constants.
4. is_constant(A) trait for whether the operator is constant or not.
5. Optional: diagonal, symmetric, etc traits from LinearMaps.jl.
6. Optional: exp(A). Required for simple exponential integration.
7. Optional: expmv(A,u,p,t) = exp(t*A)*u and expmv!(v,A::DiffEqOperator,u,p,t)
   Required for sparse-saving exponential integration.
8. Optional: factorizations. A_ldiv_B, factorize et. al. This is only required
   for algorithms which use the factorization of the operator (Crank-Nicholson),
   and only for when the default linear solve is used.
=#

# Extra standard assumptions
is_constant(L::AbstractDiffEqLinearOperator) = true
# Other ones from LinearMaps.jl
# Generic fallbacks
LinearAlgebra.exp(L::AbstractDiffEqLinearOperator,t) = exp(t*L)
has_exp(L::AbstractDiffEqLinearOperator) = true
expmv(L::AbstractDiffEqLinearOperator,u,p,t) = exp(L,t)*u
expmv!(v,L::AbstractDiffEqLinearOperator,u,p,t) = A_mul_B!(v,exp(L,t),u)
# Factorizations have no fallback and just error

"""
AffineDiffEqOperator{T} <: AbstractDiffEqOperator{T}

`Ex: (A₁(t) + ... + Aₙ(t))*u + B₁(t) + ... + Bₙ(t)`

AffineDiffEqOperator{T}(As,Bs,u_cache=nothing)

Takes in two tuples for split Affine DiffEqs

1. update_coefficients! works by updating the coefficients of the component
   operators.
2. Function calls L(u,p,t) and L(u,p,t,du) are fallbacks interpretted in this form.
   This will allow them to work directly in the nonlinear ODE solvers without
   modification.
3. f(u,p,t,du) is only allowed if a u_cache is given
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
        all([size(a) == size(As[1])
             for a in As]) || error("Operator sizes do not agree")
        new{T,typeof(As),typeof(Bs),typeof(u_cache)}(As,Bs,u_cache)
    end
end

Base.size(L::AffineDiffEqOperator) = size(L.As[1])


function (L::AffineDiffEqOperator)(u,p,t::Number)
    tmp = sum((update_coefficients(A,u,p,t); A*u for A in L.As))
    tmp2 = sum((typeof(B) <: Union{Number,AbstractArray} ? B : B(t) for B in L.Bs))
    tmp + tmp2
end

function (L::AffineDiffEqOperator)(du,u,p,t::Number)
    update_coefficients!(L,u,p,t)
    L.u_cache == nothing && error("Can only use inplace AffineDiffEqOperator if u_cache is given.")
    u_cache = L.u_cache
    fill!(du,zero(first(du)))
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

function update_coefficients!(L::AffineDiffEqOperator,u,p,t)
    # TODO: Make type-stable via recursion
    for A in L.As; update_coefficients!(A,u,p,t); end
    for B in L.Bs; update_coefficients!(B,u,p,t); end
end
