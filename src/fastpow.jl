# From David Goldberg's blog post
# https://tech.ebayinc.com/engineering/fast-approximate-logarithms-part-iii-the-formulas/
@inline function fastlog2(x::Float32)::Float32
    # (x-1)*(a*(x-1) + b)/((x-1) + c) (line 8 of table 2)
    a =   0.338953f0
    b =   2.198599f0
    c =   1.523692f0
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
        fexp = exp - 126f0    # 126 instead of 127 compensates for division by 2
        signif = signif - 1.0f0
    else
        ux2i = (ux1i & 0x007FFFFF) | 0x3f800000
        signif = reinterpret(Float32, ux2i)
        fexp = exp - 127f0
        signif = signif - 1.0f0
    end
    lg2 = fexp + signif*(a*signif + b)/(signif + c)
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
@inline function _exp2(x::Float32)
    TBLBITS = UInt32(4)
    TBLSIZE = UInt32(1 << TBLBITS)

    redux = Float32(0x1.8p23f) / TBLSIZE
    P1    = Float32(0x1.62e430p-1f)
    P2    = Float32(0x1.ebfbe0p-3f)
    P3    = Float32(0x1.c6b348p-5f)
    P4    = Float32(0x1.3b2c9cp-7f)

    # Reduce x, computing z, i0, and k.
    t::Float32 = x + redux
    i0 = reinterpret(UInt32, t)
    i0 += TBLSIZE ÷ UInt32(2)
    k::UInt32 = unsafe_trunc(UInt32, (i0 >> TBLBITS) << 20)
    i0 &= TBLSIZE - UInt32(1)
    t -= redux
    z = x - t
    twopk = Float32(reinterpret(Float64, UInt64(0x3ff00000 + k) << 32))

    # Compute r = exp2(y) = exp2ft[i0] * p(z).
    tv = EXP2FT[i0+UInt32(1)]
    u = tv * z
    tv = tv + u * (P1 + z * P2) + u * (z * z) * (P3 + z * P4)

    # Scale by 2**(k>>20)
    return tv * twopk
end

"""
    fastpow(x::Real, y::Real) -> Float32
"""
@inline function fastpow(x::Real, y::Real)
    if iszero(x)
        return 0f0
    elseif isinf(x) && isinf(y)
        return Float32(Inf)
    else
        return _exp2(convert(Float32,y) * fastlog2(convert(Float32, x)))
    end
end
@inline fastpow(x, y) = x^y
