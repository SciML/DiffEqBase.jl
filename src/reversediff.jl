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
