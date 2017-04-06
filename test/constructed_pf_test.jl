using DiffEqBase
using Base.Test

pf_func = function (t,u,p,du)
  du[1] = p[1] * u[1] - p[2] * u[1]*u[2]
  du[2] = -3 * u[2] + u[1]*u[2]
end

pf = ParameterizedFunction(pf_func,[1.5,1.0])

pf_func2 = function (t,u,p)
  [p[1] * u[1] - p[2] * u[1]*u[2];-3 * u[2] + u[1]*u[2]]
end

pf2 = ParameterizedFunction(pf_func2,[1.5,1.0])

@test param_values(pf) == Real[1.5,1]
@test num_params(pf) == 2

t = 1.0
u = [2.0,3.0]
du = zeros(2)
pf(t,u,du)

pf(t,u,du)
@test du == [-3.0,-3.0]
@test pf2(t,u) == [-3.0,-3.0]
