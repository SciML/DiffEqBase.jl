ZygoteRules.@adjoint ODESolution(u,args...) = ODESolution(u,args...), ȳ -> (ȳ,ntuple(_->nothing, length(args))...)

ZygoteRules.@adjoint function ODESolution{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11
                }(u,args...) where {T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11}

                f = function (ȳ)
                  (ȳ,ntuple(_->nothing, length(args))...)
                end

                ODESolution{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11}(u,args...),f
end

ZygoteRules.@adjoint function getindex(sol::DESolution, i)
  sol[i], function (Δ)
    Δ′ = Union{Nothing, eltype(sol.u)}[nothing for x in sol.u]
    Δ′[i] = Δ
    (Δ′,nothing)
  end
end

ZygoteRules.@adjoint function getindex(sol::DESolution, i, j...)
  sol[i,j...], function (Δ)
    Δ′ = zero(sol)
    Δ′[i,j...] = Δ
    (Δ′, map(_ -> nothing, i)...)
  end
end
