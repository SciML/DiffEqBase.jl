"""
AffineDiffEqOperator{T} <: AbstractDiffEqOperator{T}

`Ex: (A₁(t) + ... + Aₙ(t))*u + B₁(t) + ... + Bₙ(t)`

AffineDiffEqOperator{T}(As,Bs,du_cache=nothing)

Takes in two tuples for split Affine DiffEqs

1. update_coefficients! works by updating the coefficients of the component
   operators.
2. Function calls L(u, p, t) and L(du, u, p, t) are fallbacks interpretted in this form.
   This will allow them to work directly in the nonlinear ODE solvers without
   modification.
3. f(du, u, p, t) is only allowed if a du_cache is given
4. B(t) can be Union{Number,AbstractArray}, in which case they are constants.
   Otherwise they are interpreted they are functions v=B(t) and B(v,t)

Solvers will see this operator from integrator.f and can interpret it by
checking the internals of As and Bs. For example, it can check isconstant(As[1])
etc.
"""
struct AffineDiffEqOperator{T,T1,T2,U} <: AbstractDiffEqOperator{T}
    As::T1
    Bs::T2
    du_cache::U
    function AffineDiffEqOperator{T}(As,Bs,du_cache=nothing) where T
        all([size(a) == size(As[1])
             for a in As]) || error("Operator sizes do not agree")
        new{T,typeof(As),typeof(Bs),typeof(du_cache)}(As,Bs,du_cache)
    end
end

Base.size(L::AffineDiffEqOperator) = size(L.As[1])


function (L::AffineDiffEqOperator)(u,p,t::Number)
    update_coefficients!(L,u,p,t)
    du = sum(A*u for A in L.As)
    for B in L.Bs
        if typeof(B) <: Union{Number,AbstractArray}
            du .+= B
        else
            du .+= B(t)
        end
    end
    du
end

function (L::AffineDiffEqOperator)(du,u,p,t::Number)
    update_coefficients!(L,u,p,t)
    L.du_cache === nothing && error("Can only use inplace AffineDiffEqOperator if du_cache is given.")
    du_cache = L.du_cache
    fill!(du,zero(first(du)))
    # TODO: Make type-stable via recursion
    for A in L.As
        mul!(du_cache,A,u)
        du .+= du_cache
    end
    for B in L.Bs
        if typeof(B) <: Union{Number,AbstractArray}
            du .+= B
        else
            B(du_cache,t)
            du .+= du_cache
        end
    end
end

function update_coefficients!(L::AffineDiffEqOperator,u,p,t)
    # TODO: Make type-stable via recursion
    for A in L.As; update_coefficients!(A,u,p,t); end
    for B in L.Bs; update_coefficients!(B,u,p,t); end
end

@deprecate is_constant(L::AbstractDiffEqOperator) isconstant(L)
