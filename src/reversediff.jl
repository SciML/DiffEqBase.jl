function concrete_solve(prob::DiffEqBase.DEProblem,alg::Union{DiffEqBase.DEAlgorithm,Nothing},u0::ReverseDiff.TrackedArray,p::ReverseDiff.TrackedArray,args...;
                                      sensealg=nothing,kwargs...)
  ReverseDiff.track(concrete_solve,prob,alg,u0,p,args...;sensealg=sensealg,kwargs...)
end

function concrete_solve(prob::DiffEqBase.DEProblem,alg::Union{DiffEqBase.DEAlgorithm,Nothing},u0,p::ReverseDiff.TrackedArray,args...;
                                      sensealg=nothing,kwargs...)
  ReverseDiff.track(concrete_solve,prob,alg,u0,p,args...;sensealg=sensealg,kwargs...)
end

function concrete_solve(prob::DiffEqBase.DEProblem,alg::Union{DiffEqBase.DEAlgorithm,Nothing},u0::ReverseDiff.TrackedArray,p,args...;
                                      sensealg=nothing,kwargs...)
  ReverseDiff.track(concrete_solve,prob,alg,u0,p,args...;sensealg=sensealg,kwargs...)
end

ReverseDiff.@grad function concrete_solve(prob,alg,u0,p,args...;
                                          sensealg=nothing,kwargs...)
  out = _concrete_solve_adjoint(prob,alg,sensealg,ReverseDiff.value(u0),ReverseDiff.value(p),args...;kwargs...)
  Array(out[1]),out[2]
end
