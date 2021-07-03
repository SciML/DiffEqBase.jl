using LinearAlgebra, SparseArrays, SuiteSparse

struct ForwardSensitivityW{T,TW1<:AbstractMatrix{T},TW2} <: AbstractMatrix{T}
  W1::TW1
  W2::TW2
end

struct ForwardSensitivityWFactorization{T,TF1<:Factorization{T},TF2} <: Factorization{T}
  F1::TF1
  F2::TF2
end
for ff in [:lu!, :qr!, :cholesky!, :svd!], f in [ff, Symbol(string(ff)[1:end-1])]
  @eval function LinearAlgebra.$f(W::ForwardSensitivityW)
    fact1 = $f(W.W1)
    fact2 = W.W2 === nothing ? W.W2 : $f(W.W2)
    ForwardSensitivityWFactorization(fact1, fact2)
  end
end
function LinearAlgebra.ldiv!(F::ForwardSensitivityWFactorization, x)
  @unpack F1, F2 = F
  n = size(F1, 1)
  k = length(x)÷n
  @assert k*n == length(x) && k > 1
  if F2 === nothing
    ldiv!(F1, reshape(x, n, k))
  else
    ldiv!(F1, @view x[1:n])
    ldiv!(F2, reshape((@view x[n+1:end]), n, k-1))
  end
  x
end
function Base.:(\)(W::ForwardSensitivityW, x)
  @unpack W1, W2 = W
  n = size(W1, 1)
  k = length(x)÷n
  @assert k*n == length(x) && k > 1
  if W2 === nothing
    return W \ reshape(x, n, k)
  else
    [
     F1 \ (@view x[1:n])
     F2 \ reshape((@view x[n+1:end]), n, k-1)
    ]
  end
end

# This is only used for oop stiff solvers
default_factorize(A) = lu(A)

mutable struct LinSolveFactorize{F}
  factorization::F
  A
end
LinSolveFactorize(factorization) = LinSolveFactorize(factorization,nothing)
function (p::LinSolveFactorize)(x,A,b,update_matrix=false;kwargs...)
  if update_matrix
    p.A = p.factorization(A)
  end
  if typeof(p.A) <: SuiteSparse.UMFPACK.UmfpackLU || typeof(p.factorization) <: typeof(lu)
    ldiv!(x,p.A,b) # No 2-arg form for SparseArrays!
  else
    x .= b
    ldiv!(p.A,x)
  end
end
function (p::LinSolveFactorize)(::Type{Val{:init}},f,u0_prototype)
  LinSolveFactorize(p.factorization,nothing)
end

mutable struct LinSolveGPUFactorize{F,T}
  factorization::F
  A
  x_cache::T
end
LinSolveGPUFactorize(factorization=qr) = LinSolveGPUFactorize(factorization,nothing,nothing)
function (p::LinSolveGPUFactorize)(x,A,b,update_matrix=false;kwargs...)
  if update_matrix
    p.A = p.factorization(cuify(A))
  end
  ldiv!(p.x_cache,p.A,cuify(b))
  x .= Array(p.x_cache)
end
function (p::LinSolveGPUFactorize)(::Type{Val{:init}},f,u0_prototype)
  LinSolveGPUFactorize(p.factorization,nothing,cuify(u0_prototype))
end

### Default Linsolve

# Try to be as smart as possible
# lu! if Matrix
# lu if sparse
# gmres if operator

mutable struct DefaultLinSolve
  A
  iterable
end
DefaultLinSolve() = DefaultLinSolve(nothing, nothing)
@noinline function checkreltol(reltol)
  if !(reltol isa Real)
    error("Non real valued reltol is not supported by the linear iterative solvers. To customize tolerances for the linear iterative solvers, use the syntax like `KenCarp3(linsolve=LinSolveGMRES(abstol=1e-16,reltol=1e-16))`.")
  end
  return reltol
end

function (p::DefaultLinSolve)(x,A,b,update_matrix=false;reltol=nothing, kwargs...)
  if p.iterable isa Vector && eltype(p.iterable) <: LinearAlgebra.BlasInt # `iterable` here is the pivoting vector
    F = LU{eltype(A)}(A, p.iterable, zero(LinearAlgebra.BlasInt))
    ldiv!(x, F, b)
    return nothing
  end
  if update_matrix
    if A isa Matrix
      blasvendor = BLAS.vendor()
      # if the user doesn't use OpenBLAS, we assume that is a better BLAS
      # implementation like MKL
      #
      # RecursiveFactorization seems to be consistantly winning below 100
      # https://discourse.julialang.org/t/ann-recursivefactorization-jl/39213
      if ArrayInterface.can_setindex(x) && (size(A,1) <= 100 || ((blasvendor === :openblas || blasvendor === :openblas64) && size(A,1) <= 500))
        p.A = RecursiveFactorization.lu!(A)
      else
        p.A = lu!(A)
      end
    elseif A isa Union{Tridiagonal, ForwardSensitivityW}
      p.A = lu!(A)
    elseif A isa Union{SymTridiagonal}
      p.A = ldlt!(A)
    elseif A isa Union{Symmetric,Hermitian}
      p.A = bunchkaufman!(A)
    elseif A isa SparseMatrixCSC
      p.A = lu(A)
    elseif ArrayInterface.isstructured(A)
      p.A = factorize(A)
    elseif !(A isa AbstractDiffEqOperator)
      # Most likely QR is the one that is overloaded
      # Works on things like CuArrays
      p.A = qr(A)
    end
  end

  if A isa Union{Matrix,SymTridiagonal,Tridiagonal,Symmetric,Hermitian,ForwardSensitivityW} # No 2-arg form for SparseArrays!
    x .= b
    ldiv!(p.A,x)
  # Missing a little bit of efficiency in a rare case
  #elseif A isa DiffEqArrayOperator
  #  ldiv!(x,p.A,b)
  elseif ArrayInterface.isstructured(A) || A isa SparseMatrixCSC
    ldiv!(x,p.A,b)
  elseif A isa AbstractDiffEqOperator
    # No good starting guess, so guess zero
    if p.iterable === nothing
      reltol = checkreltol(reltol)
      p.iterable = IterativeSolvers.gmres_iterable!(x,A,b;initially_zero=true,restart=5,maxiter=5,abstol=1e-16,reltol=reltol,kwargs...)
    end
    x .= false
    iter = p.iterable
    purge_history!(iter, x, b)

    for residual in iter
    end
  else
    ldiv!(x,p.A,b)
  end
  return nothing
end

function (p::DefaultLinSolve)(::Type{Val{:init}},f,u0_prototype)
  if has_Wfact(f) || has_Wfact_t(f)
    piv = collect(one(LinearAlgebra.BlasInt):convert(LinearAlgebra.BlasInt, length(u0_prototype))) # pivoting vector
    DefaultLinSolve(f, piv)
  else
    DefaultLinSolve()
  end
end

const DEFAULT_LINSOLVE = DefaultLinSolve()

### Default GMRES

# Easily change to GMRES

mutable struct LinSolveIterativeSolvers{F,PL,PR,AR,A}
  generate_iterator::F
  iterable
  Pl::PL
  Pr::PR
  args::AR
  kwargs::A
end
LinSolveIterativeSolvers(generate_iterator,args...;
                         Pl=IterativeSolvers.Identity(),
                         Pr=IterativeSolvers.Identity(),
                         kwargs...) =
                         LinSolveIterativeSolvers(generate_iterator, nothing, Pl, Pr, args, kwargs)

LinSolveGMRES(args...;kwargs...) = LinSolveIterativeSolvers(IterativeSolvers.gmres_iterable!,args...;kwargs...)
LinSolveCG(args...;kwargs...) = LinSolveIterativeSolvers(IterativeSolvers.cg_iterator!,args...;kwargs...)
LinSolveBiCGStabl(args...;kwargs...) = LinSolveIterativeSolvers(IterativeSolvers.bicgstabl_iterator!,args...;kwargs...)
LinSolveChebyshev(args...;kwargs...) = LinSolveIterativeSolvers(IterativeSolvers.chebyshev_iterable!,args...;kwargs...)
LinSolveMINRES(args...;kwargs...) = LinSolveIterativeSolvers(IterativeSolvers.minres_iterable!,args...;kwargs...)

function (f::LinSolveIterativeSolvers)(x,A,b,update_matrix=false; Pl=nothing, Pr=nothing, reltol=eps(eltype(x)), kwargs...)
  if f.iterable === nothing
    Pl = ComposePreconditioner(f.Pl, Pl, true)
    Pr = ComposePreconditioner(f.Pr, Pr, false)

    reltol = checkreltol(reltol)
    f.iterable = f.generate_iterator(x,A,b,f.args...;
                                     initially_zero=true,restart=5,
                                     maxiter=5,abstol=1e-16,reltol=reltol,
                                     Pl=Pl,Pr=Pr,
                                     kwargs...,f.kwargs...)
  end
  x .= false
  iter = f.iterable
  purge_history!(iter, x, b)

  for residual in iter
  end
  return nothing
end

function (p::LinSolveIterativeSolvers)(::Type{Val{:init}},f,u0_prototype)
  LinSolveIterativeSolvers(p.generate_iterator, nothing, p.Pl, p.Pr, p.args, p.kwargs)
end

# scaling for iterative solvers
struct ScaleVector{T}
  x::T
  isleft::Bool
end
function LinearAlgebra.ldiv!(v::ScaleVector, x)
  vx = vec(v.x)
  if v.isleft
    return @.. x = x * vx
  else
    return @.. x = x / vx
  end
end
function LinearAlgebra.ldiv!(y, v::ScaleVector, x)
  vx = vec(v.x)
  if v.isleft
    return @.. y = x * vx
  else
    return @.. y = x / vx
  end
end

struct ComposePreconditioner{P, S<:Union{ScaleVector,Nothing}}
  P::P
  scale::S
  isleft::Bool
end

function LinearAlgebra.ldiv!(v::ComposePreconditioner, x)
  isid = v.P isa IterativeSolvers.Identity
  isid || ldiv!(v.P, x)
  s = v.scale
  s === nothing && return x
  ldiv!(s, x)
  return x
end

function LinearAlgebra.ldiv!(y, v::ComposePreconditioner, x)
  isid = v.P isa IterativeSolvers.Identity
  isid || ldiv!(y, v.P, x)
  s = v.scale
  s === nothing && return x
  if isid
    ldiv!(y, s, x)
  else
    if v.isleft
      @.. y = y * s.x
    else
      @.. y = y / s.x
    end
  end
  return y
end

function purge_history!(iter::IterativeSolvers.GMRESIterable, x, b)
  iter.k = 1
  iter.x  = x
  iter.b  = b

  iter.residual.current = IterativeSolvers.init!(iter.arnoldi, iter.x, iter.b, iter.Pl, iter.Ax, initially_zero = true)
  IterativeSolvers.init_residual!(iter.residual, iter.residual.current)
  iter.β = iter.residual.current
  nothing
end
