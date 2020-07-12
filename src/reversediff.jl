value(x::ReverseDiff.TrackedReal)  = x.value
value(x::ReverseDiff.TrackedArray) = x.value

# Support adaptive with non-tracked time
@inline function ODE_DEFAULT_NORM(u::ReverseDiff.TrackedArray,t) where {N}
  sqrt(sum(abs2,value(u)) / length(u))
end
@inline function ODE_DEFAULT_NORM(u::AbstractArray{<:ReverseDiff.TrackedReal,N},t) where {N}
  sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(t))) / length(u))
end
@inline function ODE_DEFAULT_NORM(u::Array{<:ReverseDiff.TrackedReal,N},t) where {N}
  sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip((value(x) for x in u),Iterators.repeated(t))) / length(u))
end
@inline ODE_DEFAULT_NORM(u::ReverseDiff.TrackedReal,t) = abs(value(u))

# Support TrackedReal time, don't drop tracking on the adaptivity there
@inline function ODE_DEFAULT_NORM(u::ReverseDiff.TrackedArray,t::ReverseDiff.TrackedReal) where {N}
  sqrt(sum(abs2,u) / length(u))
end
@inline function ODE_DEFAULT_NORM(u::AbstractArray{<:ReverseDiff.TrackedReal,N},t::ReverseDiff.TrackedReal) where {N}
  sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip(u,Iterators.repeated(t))) / length(u))
end
@inline function ODE_DEFAULT_NORM(u::Array{<:ReverseDiff.TrackedReal,N},t::ReverseDiff.TrackedReal) where {N}
  sqrt(sum(x->ODE_DEFAULT_NORM(x[1],x[2]),zip(u,Iterators.repeated(t))) / length(u))
end
@inline ODE_DEFAULT_NORM(u::ReverseDiff.TrackedReal,t::ReverseDiff.TrackedReal) = abs(u)

function solve_up(prob::DiffEqBase.DEProblem,sensealg::Union{AbstractSensitivityAlgorithm,Nothing},u0::ReverseDiff.TrackedArray,p::ReverseDiff.TrackedArray,args...;kwargs...)
  ReverseDiff.track(solve_up,prob,sensealg,u0,p,args...;kwargs...)
end

function solve_up(prob::DiffEqBase.DEProblem,sensealg::Union{AbstractSensitivityAlgorithm,Nothing},u0,p::ReverseDiff.TrackedArray,args...;kwargs...)
  ReverseDiff.track(solve_up,prob,sensealg,u0,p,args...;kwargs...)
end

function solve_up(prob::DiffEqBase.DEProblem,sensealg::Union{AbstractSensitivityAlgorithm,Nothing},u0::ReverseDiff.TrackedArray,p,args...;kwargs...)
  ReverseDiff.track(solve_up,prob,sensealg,u0,p,args...;kwargs...)
end

ReverseDiff.@grad function solve_up(prob,sensealg,u0,p,args...;kwargs...)
  out = _solve_adjoint(prob,sensealg,ReverseDiff.value(u0),ReverseDiff.value(p),args...;kwargs...)
  Array(out[1]),out[2]
end
