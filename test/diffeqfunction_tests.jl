using DiffEqBase

f(t,u,du) = (du.=u)
DE_f = DiffEqFunction{true}(f)
t = 0.0; u = Float64[1,2,3]; du = zeros(3)
DE_f(t,u,du)

NS_f = DiffEqBase.NSODEFunction(f,t,u)
NS_f = DiffEqBase.NSODEFunction{true}(f,t,u)
NS_f(t,u,du)
