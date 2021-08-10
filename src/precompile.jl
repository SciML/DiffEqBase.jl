let
    while true
        _testf(du,u,p,t) = copyto!(du,u)
        b = rand(1); x = rand(1)
        _linsolve = DEFAULT_LINSOLVE(Val{:init},ODEFunction(_testf),b)
        A = rand(1,1)
        _linsolve(x,A,b,true)
        _linsolve(x,A,b,false)
        break
    end
end
