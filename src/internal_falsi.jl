"""
  prevfloat_tdir(x, x0, x1)

Move `x` one floating point towards x0.
"""
function prevfloat_tdir(x, x0, x1)
    x1 > x0 ? prevfloat(x) : nextfloat(x)
end

function nextfloat_tdir(x, x0, x1)
    x1 > x0 ? nextfloat(x) : prevfloat(x)
end

function max_tdir(a, b, x0, x1)
    x1 > x0 ? max(a, b) : min(a, b)
end

"""
`InternalFalsi`: A non-allocating regula falsi method, internal to DiffEqBase for
simpler dependencies.
"""
struct InternalFalsi end

function SciMLBase.solve(prob::IntervalNonlinearProblem, alg::InternalFalsi, args...;
    maxiters = 1000,
    kwargs...)
    f = Base.Fix2(prob.f, prob.p)
    left, right = prob.tspan
    fl, fr = f(left), f(right)

    if iszero(fl)
        return SciMLBase.build_solution(prob, alg, left, fl;
            retcode = ReturnCode.ExactSolutionLeft, left = left,
            right = right)
    end

    i = 1
    if !iszero(fr)
        using_falsi_steps = true
        while i < maxiters
            # First, perform a regula falsi iteration
            if using_falsi_steps
                if nextfloat_tdir(left, prob.tspan...) == right
                    return SciMLBase.build_solution(prob, alg, left, fl;
                        retcode = ReturnCode.FloatingPointLimit,
                        left = left, right = right)
                end
                mid = (fr * left - fl * right) / (fr - fl)
                for i in 1:10
                    mid = max_tdir(left, prevfloat_tdir(mid, prob.tspan...), prob.tspan...)
                end
                if mid == right || mid == left
                    using_falsi_steps = false
                    continue
                end
                fm = f(mid)
                if iszero(fm)
                    right = mid
                    using_falsi_steps = false
                    continue
                end
                if sign(fl) == sign(fm)
                    fl = fm
                    left = mid
                else
                    fr = fm
                    right = mid
                end
                i += 1
            end

            # Then, perform a bisection iteration
            mid = (left + right) / 2
            (mid == left || mid == right) &&
                return SciMLBase.build_solution(prob, alg, left, fl;
                    retcode = ReturnCode.FloatingPointLimit,
                    left = left, right = right)
            fm = f(mid)
            if iszero(fm)
                right = mid
                fr = fm
            elseif sign(fm) == sign(fl)
                left = mid
                fl = fm
            else
                right = mid
                fr = fm
            end
            i += 1
        end
    end

    return SciMLBase.build_solution(prob, alg, left, fl; retcode = ReturnCode.MaxIters,
        left = left, right = right)
end

function scalar_nlsolve_ad(prob, alg::InternalFalsi, args...; kwargs...)
    f = prob.f
    p = value(prob.p)

    if prob isa IntervalNonlinearProblem
        tspan = value(prob.tspan)
        newprob = IntervalNonlinearProblem(f, tspan, p; prob.kwargs...)
    else
        u0 = value(prob.u0)
        newprob = NonlinearProblem(f, u0, p; prob.kwargs...)
    end

    sol = solve(newprob, alg, args...; kwargs...)

    uu = sol.u
    if p isa Number
        f_p = ForwardDiff.derivative(Base.Fix1(f, uu), p)
    else
        f_p = ForwardDiff.gradient(Base.Fix1(f, uu), p)
    end

    f_x = ForwardDiff.derivative(Base.Fix2(f, p), uu)
    pp = prob.p
    sumfun = let f_x′ = -f_x
        ((fp, p),) -> (fp / f_x′) * ForwardDiff.partials(p)
    end
    partials = sum(sumfun, zip(f_p, pp))
    return sol, partials
end

function SciMLBase.solve(prob::IntervalNonlinearProblem{uType, iip,
    <:ForwardDiff.Dual{T, V, P}},
    alg::InternalFalsi, args...;
    kwargs...) where {uType, iip, T, V, P}
    sol, partials = scalar_nlsolve_ad(prob, alg, args...; kwargs...)
    return SciMLBase.build_solution(prob, alg, Dual{T, V, P}(sol.u, partials),
        sol.resid; retcode = sol.retcode,
        left = Dual{T, V, P}(sol.left, partials),
        right = Dual{T, V, P}(sol.right, partials))
end

function SciMLBase.solve(prob::IntervalNonlinearProblem{uType, iip,
    <:AbstractArray{
        <:ForwardDiff.Dual{T,
            V,
            P},
    }},
    alg::InternalFalsi, args...;
    kwargs...) where {uType, iip, T, V, P}
    sol, partials = scalar_nlsolve_ad(prob, alg, args...; kwargs...)
    
    return SciMLBase.build_solution(prob, alg, Dual{T, V, P}(sol.u, partials),
        sol.resid; retcode = sol.retcode,
        left = Dual{T, V, P}(sol.left, partials),
        right = Dual{T, V, P}(sol.right, partials))
end
