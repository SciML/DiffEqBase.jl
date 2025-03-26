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
`InternalITP`: A non-allocating ITP method, internal to DiffEqBase for
simpler dependencies.
"""
struct InternalITP
    scaled_k1::Float64
    n0::Int
end

InternalITP() = InternalITP(0.2, 10)

function SciMLBase.solve(prob::IntervalNonlinearProblem{IP, Tuple{T, T}}, alg::InternalITP,
        args...;
        maxiters = 1000, kwargs...) where {IP, T}
    f = Base.Fix2(prob.f, prob.p)
    left, right = prob.tspan # a and b
    fl, fr = f(left), f(right)
    ϵ = eps(T)
    if iszero(fl)
        return SciMLBase.build_solution(prob, alg, left, fl;
            retcode = ReturnCode.ExactSolutionLeft, left, right)
    elseif iszero(fr)
        return SciMLBase.build_solution(prob, alg, right, fr;
            retcode = ReturnCode.ExactSolutionRight, left, right)
    end
    span = abs(right - left)
    k1 = T(alg.scaled_k1)/span
    n0 = T(alg.n0)
    n_h = exponent(span / (2 * ϵ))
    ϵ_s = ϵ * exp2(n_h + n0)
    T0 = zero(fl)

    i = 1
    while i ≤ maxiters
        stats.nsteps += 1
        span = abs(right - left)
        mid = (left + right) / 2
        r = ϵ_s - (span / 2)

        x_f = left + span * (fl / (fl - fr))  # Interpolation Step

        δ = max(k1 * span^2, eps(x_f))
        diff = mid - x_f

        xt = ifelse(δ ≤ abs(diff), x_f + copysign(δ, diff), mid)  # Truncation Step

        xp = ifelse(abs(xt - mid) ≤ r, xt, mid - copysign(r, diff))  # Projection Step
        if span < 2ϵ
            return SciMLBase.build_solution(
                prob, alg, xt, f(xt); retcode = ReturnCode.Success, left, right
            )
        end
        yp = f(xp)
        stats.nf += 1
        yps = yp * sign(fr)
        if yps > T0
            right, fr = xp, yp
        elseif yps < T0
            left, fl = xp, yp
        else
            return SciMLBase.build_solution(
                prob, alg, xp, yps; retcode = ReturnCode.Success, left, right
            )
        end

        i += 1
        ϵ_s /= 2

        if nextfloat_tdir(left, left, right) == right
            return SciMLBase.build_solution(
                prob, alg, right, fr; retcode = ReturnCode.FloatingPointLimit, left, right
            )
        end
    end
    return SciMLBase.build_solution(prob, alg, left, fl; retcode = ReturnCode.MaxIters,
        left = left, right = right)
end
