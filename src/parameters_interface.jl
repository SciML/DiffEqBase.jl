function param_values(f::AbstractParameterizedFunction)
  [getfield(f,s) for s in f.params]
end

num_params(f::AbstractParameterizedFunction) = length(f.params)

# Fallbacks

param_values(f) = nothing
num_params(f) = 0

function problem_new_parameters(prob::ODEProblem,p)
  f = (t,u,du) -> prob.f(t,u,p,du)
  uEltype = eltype(p)
  u0 = [uEltype(prob.u0[i]) for i in 1:length(prob.u0)]
  tspan = (uEltype(prob.tspan[1]),uEltype(prob.tspan[2]))
  ODEProblem(f,u0,tspan)
end
param_values(prob::ODEProblem) = param_values(prob.f)
num_params(prob::ODEProblem) = num_params(prob.f)

function problem_new_parameters(prob::SDEProblem,p)
  fpars = num_params(prob.f)
  f = (t,u,du) -> prob.f(t,u,@view(p[1:fpars]),du)
  if num_params(prob.g) != 0
    g = (t,u,du) -> prob.g(t,u,@view(p[(fpars+1):end],du)
  else
    g = prob.g
  end
  uEltype = eltype(p)
  u0 = [uEltype(prob.u0[i]) for i in 1:length(prob.u0)]
  tspan = (uEltype(prob.tspan[1]),uEltype(prob.tspan[2]))
  ODEProblem(f,g,u0,tspan)
end
param_values(prob::SDEProblem) = [param_values(prob.f) ; param_values(prob.g)]
num_params(prob::SDEProblem) = num_params(prob.f) + num_params(prob.g)
