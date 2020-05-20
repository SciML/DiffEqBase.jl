value(x::Type{Tracker.TrackedReal{T}}) where T = T
value(x::Type{Tracker.TrackedArray{T,N,A}}) where {T,N,A} = Array{T,N}
value(x::Tracker.TrackedReal)  = x.data
value(x::Tracker.TrackedArray) = x.data

@inline fastpow(x::Tracker.TrackedReal, y::Tracker.TrackedReal) = x^y
@inline Base.any(f::Function,x::Tracker.TrackedArray) = any(f,Tracker.data(x))

# Support adaptive with non-tracked time
@inline function ODE_DEFAULT_NORM(u::Tracker.TrackedArray,t) where {N}
  sqrt(sum(abs2,value(u)) / length(u))
end
@inline function ODE_DEFAULT_NORM(u::AbstractArray{<:Tracker.TrackedReal,N},t) where {N}
  sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(t))) / length(u))
end
@inline function ODE_DEFAULT_NORM(u::Array{<:Tracker.TrackedReal,N},t) where {N}
  sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(t))) / length(u))
end
@inline ODE_DEFAULT_NORM(u::Tracker.TrackedReal,t) = abs(value(u))

# Support TrackedReal time, don't drop tracking on the adaptivity there
@inline function ODE_DEFAULT_NORM(u::Tracker.TrackedArray,t::Tracker.TrackedReal) where {N}
  sqrt(sum(abs2,u) / length(u))
end
@inline function ODE_DEFAULT_NORM(u::AbstractArray{<:Tracker.TrackedReal,N},t::Tracker.TrackedReal) where {N}
  sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip(u,Iterators.repeated(t))) / length(u))
end
@inline function ODE_DEFAULT_NORM(u::Array{<:Tracker.TrackedReal,N},t::Tracker.TrackedReal) where {N}
  sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip(u,Iterators.repeated(t))) / length(u))
end
@inline ODE_DEFAULT_NORM(u::Tracker.TrackedReal,t::Tracker.TrackedReal) = abs(u)

function DiffEqBase.concrete_solve(prob::DiffEqBase.DEProblem,alg::DiffEqBase.DEAlgorithm,u0::Tracker.TrackedArray,p::Tracker.TrackedArray,args...;
                                      sensealg=nothing,kwargs...)
  Tracker.track(concrete_solve,prob,alg,u0,p,args...;sensealg=sensealg,kwargs...)
end

function DiffEqBase.concrete_solve(prob::DiffEqBase.DEProblem,alg::DiffEqBase.DEAlgorithm,u0,p::Tracker.TrackedArray,args...;
                                      sensealg=nothing,kwargs...)
  Tracker.track(concrete_solve,prob,alg,u0,p,args...;sensealg=sensealg,kwargs...)
end

function DiffEqBase.concrete_solve(prob::DiffEqBase.DEProblem,alg::DiffEqBase.DEAlgorithm,u0::Tracker.TrackedArray,p,args...;
                                      sensealg=nothing,kwargs...)
  Tracker.track(concrete_solve,prob,alg,u0,p,args...;sensealg=sensealg,kwargs...)
end

Tracker.@grad function concrete_solve(prob,alg,u0,p,args...;
                                      sensealg=nothing,kwargs...)
  _concrete_solve_adjoint(prob,alg,sensealg,Tracker.data(u0),Tracker.data(p),args...;kwargs...)
end
