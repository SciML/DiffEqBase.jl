let
    while true
        _testf(du,u,p,t) = copyto!(du,u)
        b = rand(1); x = rand(1)
        _linsolve = DEFAULT_LINSOLVE(Val{:init},ODEFunction(_testf),b)
        A = rand(1,1)
        _linsolve(x,A,b,true)
        _linsolve(x,A,b,false)
        _linsolve = LUFactorize()(Val{:init},ODEFunction(_testf),b)
        _linsolve(x,A,b,true)
        _linsolve(x,A,b,false)
        Pl = ScaleVector([1.0],true)
        Pr = ScaleVector([1.0],false)
        reltol = 1.0
        _linsolve(x,A,b,true;reltol=reltol,Pl=Pl,Pr=Pr)
        _linsolve(x,A,b,false;reltol=reltol,Pl=Pl,Pr=Pr)        
        break
    end
end
