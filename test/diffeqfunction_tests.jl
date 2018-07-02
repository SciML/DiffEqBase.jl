using DiffEqBase, Test

f_op = (u,p,t) -> u
f_ip = (du,u,p,t) -> du .= u
DE_f = DiffEqFunction{false}(f_op)
DE_f_ip = DiffEqFunction{true}(f_ip)
du = zeros(3); u = [1.0,2.0,3.0]; p = nothing; t = 0.0
@test DE_f(u,p,t) == u
DE_f_ip(du,u,p,t)
@test du == u
