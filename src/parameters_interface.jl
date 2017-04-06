function param_values(f::AbstractParameterizedFunction)
  [getfield(f,s) for s in f.params]
end

num_params(f::AbstractParameterizedFunction) = length(f.params)

function problem_new_parameters(prob::ODEProblem,p)
  f = (t,u,du) -> prob.f(t,u,p,du)
  uEltype = eltype(p)
  u0 = [uEltype(prob.u0[i]) for i in 1:length(prob.u0)]
  tspan = (uEltype(prob.tspan[1]),uEltype(prob.tspan[2]))
  ODEProblem(f,u0,tspan)
end
