_linsolve = DefaultLinSolve()
_testf(du,u,p,t) = du .= u
b = rand(1); x = rand(1)
DEFAULT_LINSOLVE(Val{:init},ODEFunction(_testf),b)
A = rand(1,1)
DEFAULT_LINSOLVE(x,A,b,true)
