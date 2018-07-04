using DiffEqBase, Test

f_op = (u,p,t) -> u
f_ip = (du,u,p,t) -> du .= u
DE_f = ODEFunction{false}(f_op)
DE_f_ip = ODEFunction{true}(f_ip)
du = zeros(3); u = [1.0,2.0,3.0]; p = nothing; t = 0.0
@test DE_f(u,p,t) == u
DE_f_ip(du,u,p,t)
@test du == u

f(u,p,t) = u
jac = (u,p,t) -> 1
@inferred ODEFunction{false}(f_op,jac=jac)

f(du,u,p,t) = u
