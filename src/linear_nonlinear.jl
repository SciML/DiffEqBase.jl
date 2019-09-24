using LinearAlgebra, SparseArrays, SuiteSparse, KrylovKit

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

function (p::DefaultLinSolve)(x,A,b,update_matrix=false;tol=nothing, kwargs...)
  if p.iterable isa Vector && eltype(p.iterable) <: LinearAlgebra.BlasInt # `iterable` here is the pivoting vector
    F = LU{eltype(A)}(A, p.iterable, zero(LinearAlgebra.BlasInt))
    ldiv!(x, F, b)
    return nothing
  end
  if update_matrix
    if typeof(A) <: Matrix
      blasvendor = BLAS.vendor()
      if (blasvendor === :openblas || blasvendor === :openblas64) && size(A,1) <= 500 && ArrayInterface.can_setindex(x) # if the user doesn't use OpenBLAS, we assume that is a much better BLAS implementation like MKL
        p.A = RecursiveFactorization.lu!(A)
      else
        p.A = lu!(A)
      end
    elseif typeof(A) <: Tridiagonal
      p.A = lu!(A)
    elseif typeof(A) <: Union{SymTridiagonal}
      p.A = ldlt!(A)
    elseif typeof(A) <: Union{Symmetric,Hermitian}
      p.A = bunchkaufman!(A)
    elseif typeof(A) <: SparseMatrixCSC
      p.A = lu(A)
    elseif ArrayInterface.isstructured(A)
      p.A = factorize(A)
    elseif !(typeof(A) <: AbstractDiffEqOperator)
      # Most likely QR is the one that is overloaded
      # Works on things like CuArrays
      p.A = qr(A)
    end
  end

  if typeof(A) <: Union{Matrix,SymTridiagonal,Tridiagonal,Symmetric,Hermitian} # No 2-arg form for SparseArrays!
    x .= b
    ldiv!(p.A,x)
  # Missing a little bit of efficiency in a rare case
  #elseif typeof(A) <: DiffEqArrayOperator
  #  ldiv!(x,p.A,b)
  elseif ArrayInterface.isstructured(A) || A isa SparseMatrixCSC
    ldiv!(x,p.A,b)
  elseif typeof(A) <: AbstractDiffEqOperator
    # No good starting guess, so guess zero
    if p.iterable === nothing
      p.iterable = IterativeSolvers.gmres_iterable!(x,A,b;initially_zero=true,restart=5,maxiter=5,tol=1e-16,kwargs...)
      p.iterable.reltol = tol
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

function (f::LinSolveIterativeSolvers)(x,A,b,update_matrix=false; Pl=nothing, Pr=nothing, tol=nothing, kwargs...)
  if f.iterable === nothing
    Pl = ComposePreconditioner(f.Pl, Pl, true)
    Pr = ComposePreconditioner(f.Pr, Pr, false)
    f.iterable = f.generate_iterator(x,A,b,f.args...;
                                     initially_zero=true,restart=5,
                                     maxiter=5,tol=1e-16,Pl=Pl,Pr=Pr,
                                     f.kwargs...,kwargs...)
    tol′ = get(f.kwargs, :tol, nothing)
    tol′ !== nothing && (tol = tol′)
    f.iterable.reltol = tol
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

mutable struct LinSolveKrylovKit{F,PL,PR,AR,A}
  solver::F
  algo
  # Preconditioners, not yet in KrylovKit:
  Pl::PL
  Pr::PR
  args::AR
  kwargs::A
end

LinSolveKrylovKitGMRES(generate_iterator,args...;
                         Pl=IterativeSolvers.Identity(),
                         Pr=IterativeSolvers.Identity(),
                         kwargs...) =
                         LinSolveKrylovKit(KrylovKit.linsolve, GMRES(kwargs),
                         nothing, nothing,
                         args, kwargs)

LinSolveKrylovKitCG(generate_iterator,args...;
                        Pl=IterativeSolvers.Identity(),
                        Pr=IterativeSolvers.Identity(),
                        kwargs...) =
                        LinSolveKrylovKit(KrylovKit.linsolve, CG(kwargs),
                        nothing, nothing,
                        args, kwargs)

function (f::LinSolveKrylovKitGMRES)(x,A,b,update_matrix=false; Pl=nothing, Pr=nothing, tol=nothing, kwargs...)
  sol, info = f.generate_iterator(A,b,f.algo,f.args...;
                                  restart=5,
                                  maxiter=5,tol=1e-16,
                                  f.kwargs...,kwargs...)
  (info.converged == 0) && @warn "GMRES did not converged, $info"
  copyto!(x, sol)

  return nothing
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
