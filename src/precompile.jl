let
    while true
        _linsolve = DefaultLinSolve()
        _testf(du,u,p,t) = du .= u
        b = rand(1); x = rand(1)
        _linsolve(Val{:init},ODEFunction(_testf),b)
        A = rand(1,1)
        _linsolve(x,A,b,true)
        _linsolve(x,A,b,false)
        break
    end
end
