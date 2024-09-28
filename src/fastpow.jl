# From David Goldberg's blog post
# https://tech.ebayinc.com/engineering/fast-approximate-logarithms-part-iii-the-formulas/
@inline function fastlog2(x::Float32)::Float32
    # (x-1)*(a*(x-1) + b)/((x-1) + c) (line 8 of table 2)
    a = 0.338953f0
    b = 2.198599f0
    c = 1.523692f0
    #
    # Assume IEEE representation, which is sgn(1):exp(8):frac(23)
    # representing (1+frac)*2^(exp-127)  Call 1+frac the significand
    #

    # get exponent
    ux1i = reinterpret(UInt32, x)
    exp = (ux1i & 0x7F800000) >> 23
    # actual exponent is exp-127, will subtract 127 later

    greater = ux1i & 0x00400000  # true if signif > 1.5
    if greater !== 0x00000000
        ux2i = (ux1i & 0x007FFFFF) | 0x3f000000
        signif = reinterpret(Float32, ux2i)
        fexp = exp - 126.0f0    # 126 instead of 127 compensates for division by 2
        signif = signif - 1.0f0
    else
        ux2i = (ux1i & 0x007FFFFF) | 0x3f800000
        signif = reinterpret(Float32, ux2i)
        fexp = exp - 127.0f0
        signif = signif - 1.0f0
    end
    lg2 = fexp + signif * (a * signif + b) / (signif + c)
    return lg2
end

# Translated from OpenLibm but less accurate because I forced the tableau to be
# Float32, whereas OpenLibm uses Float64
#
# https://github.com/JuliaMath/openlibm/blob/cca41bc1abd01804afa4862bbd2c79cc9803171a/src/s_exp2f.c
const EXP2FT = (Float32(0x1.6a09e667f3bcdp-1),
    Float32(0x1.7a11473eb0187p-1),
    Float32(0x1.8ace5422aa0dbp-1),
    Float32(0x1.9c49182a3f090p-1),
    Float32(0x1.ae89f995ad3adp-1),
    Float32(0x1.c199bdd85529cp-1),
    Float32(0x1.d5818dcfba487p-1),
    Float32(0x1.ea4afa2a490dap-1),
    Float32(0x1.0000000000000p+0),
    Float32(0x1.0b5586cf9890fp+0),
    Float32(0x1.172b83c7d517bp+0),
    Float32(0x1.2387a6e756238p+0),
    Float32(0x1.306fe0a31b715p+0),
    Float32(0x1.3dea64c123422p+0),
    Float32(0x1.4bfdad5362a27p+0),
    Float32(0x1.5ab07dd485429p+0))

"""
    fastpow(x::T, y::T) where {T} -> float(T)
    Trips through Float32 for performance.
"""
@inline function fastpow(x::T, y::T) where {T<:Real}
    outT = float(T)
    if iszero(x)
        return zero(outT)
    elseif isinf(x) && isinf(y)
        return convert(outT, Inf)
    else
        return convert(
            outT, @fastmath exp2(convert(Float32, y) * fastlog2(convert(Float32, x))))
    end
end

@inline fastpow(x, y) = x^y
