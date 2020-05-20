ReverseDiff.@grad function concrete_solve(prob,alg,u0,p,args...;
                                          sensealg=nothing,kwargs...)
  _concrete_solve_adjoint(prob,alg,sensealg,u0,p,args...;kwargs...)
end
