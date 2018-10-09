using DiffEqBase, Test

f_op = (u,p,t) -> u
f_ip = (du,u,p,t) -> du .= u
DE_f = ODEFunction{false}(f_op)
DE_f_ip = ODEFunction{true}(f_ip)
DI_f = DiscreteFunction{false}(f_op)
DI_f_ip = DiscreteFunction{true}(f_ip)
du = zeros(3); u = [1.0,2.0,3.0]; p = nothing; t = 0.0
@test DE_f(u,p,t) == u
DE_f_ip(du,u,p,t)
@test du == u
@test DI_f(u,p,t) == u
DI_f_ip(du,u,p,t)
@test du == u

f(u,p,t) = u
jac = (u,p,t) -> 1
jacvec = (v,u,p,t) -> 1v
@inferred ODEFunction{false}(f_op,jac=jac,jacvec=jacvec)
@inferred DiscreteFunction{false}(f_op,analytic=f)

f(du,u,p,t) = u
