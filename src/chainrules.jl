function ChainRulesCore.rrule(::Type{ODEProblem}, args...; kwargs...)
    function ODEProblemAdjoint(ȳ)
        (NoTangent(), ȳ.f, ȳ.u0, ȳ.tspan, ȳ.p, ȳ.kwargs, ȳ.problem_type)
    end

    ODEProblem(args...; kwargs...), ODEProblemAdjoint
end

function ChainRulesCore.rrule(::Type{SDEProblem}, args...; kwargs...)
    function SDEProblemAdjoint(ȳ)
        (NoTangent(), ȳ.f, ȳ.g, ȳ.u0, ȳ.tspan, ȳ.p, ȳ.kwargs, ȳ.problem_type)
    end

    SDEProblem(args...; kwargs...), SDEProblemAdjoint
end

function ChainRulesCore.rrule(::Type{
                                     <:ODESolution{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10,
                                                   T11, T12
                                                   }}, u,
                              args...) where {T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11,
                                              T12}
    function ODESolutionAdjoint(ȳ)
        (NoTangent(), ȳ, ntuple(_ -> NoTangent(), length(args))...)
    end

    ODESolution{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12}(u, args...),
    ODESolutionAdjoint
end

ZygoteRules.@adjoint function ODESolution{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12
                                          }(u,
                                            args...) where {T1, T2, T3, T4, T5, T6, T7, T8,
                                                            T9, T10, T11, T12}
    function ODESolutionAdjoint(ȳ)
        (ȳ, ntuple(_ -> nothing, length(args))...)
    end

    ODESolution{T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12}(u, args...),
    ODESolutionAdjoint
end

function ChainRulesCore.rrule(::Type{
                                     <:ODESolution{uType, tType, isinplace, P, NP, F, G, K,
                                                   ND
                                                   }}, u,
                              args...) where {uType, tType, isinplace, P, NP, F, G, K, ND}
    function SDESolutionAdjoint(ȳ)
        (NoTangent(), ȳ, ntuple(_ -> NoTangent(), length(args))...)
    end

    SDESolution{uType, tType, isinplace, P, NP, F, G, K, ND}(u, args...), SDESolutionAdjoint
end

ZygoteRules.@adjoint function SDEProblem{uType, tType, isinplace, P, NP, F, G, K, ND}(u,
                                                                                      args...) where
    {uType, tType, isinplace, P, NP, F, G, K, ND}
    function SDESolutionAdjoint(ȳ)
        (ȳ, ntuple(_ -> nothing, length(args))...)
    end

    SDESolution{uType, tType, isinplace, P, NP, F, G, K, ND}(u, args...), SDESolutionAdjoint
end

ZygoteRules.@adjoint function NonlinearSolution{T, N, uType, R, P, A, O, uType2}(u,
                                                                                 args...) where {
                                                                                                 T,
                                                                                                 N,
                                                                                                 uType,
                                                                                                 R,
                                                                                                 P,
                                                                                                 A,
                                                                                                 O,
                                                                                                 uType2
                                                                                                 }
    function NonlinearSolutionAdjoint(ȳ)
        (ȳ, ntuple(_ -> nothing, length(args))...)
    end
    NonlinearSolution{T, N, uType, R, P, A, O, uType2}(u, args...), NonlinearSolutionAdjoint
end

ZygoteRules.@adjoint function ZygoteRules.literal_getproperty(sol::AbstractTimeseriesSolution,
                                                              ::Val{:u})
    function solu_adjoint(Δ)
        zerou = zero(sol.prob.u0)
        _Δ = @. ifelse(Δ == nothing, (zerou,), Δ)
        (DiffEqBase.build_solution(sol.prob, sol.alg, sol.t, _Δ),)
    end
    sol.u, solu_adjoint
end

ZygoteRules.@adjoint function ZygoteRules.literal_getproperty(sol::AbstractNoTimeSolution,
                                                              ::Val{:u})
    function solu_adjoint(Δ)
        zerou = zero(sol.prob.u0)
        _Δ = @. ifelse(Δ == nothing, zerou, Δ)
        (DiffEqBase.build_solution(sol.prob, sol.alg, _Δ, sol.resid),)
    end
    sol.u, solu_adjoint
end

ZygoteRules.@adjoint function ZygoteRules.literal_getproperty(sol::SciMLBase.OptimizationSolution,
                                                              ::Val{:u})
    function solu_adjoint(Δ)
        zerou = zero(sol.u)
        _Δ = @. ifelse(Δ == nothing, zerou, Δ)
        (DiffEqBase.build_solution(sol.cache, sol.alg, _Δ, sol.objective),)
    end
    sol.u, solu_adjoint
end

function ChainRulesCore.rrule(::DiffEqBase.EnsembleSolution, sim, time, converged)
    out = EnsembleSolution(sim, time, converged)
    function EnsembleSolution_adjoint(p̄::AbstractArray{T, N}) where {T, N}
        arrarr = [[p̄[ntuple(x -> Colon(), Val(N - 2))..., j, i]
                   for j in 1:size(p̄)[end - 1]] for i in 1:size(p̄)[end]]
        (NoTangent(), EnsembleSolution(arrarr, 0.0, true), NoTangent(), NoTangent())
    end
    function EnsembleSolution_adjoint(p̄::EnsembleSolution)
        (NoTangent(), p̄, NoTangent(), NoTangent())
    end
    out, EnsembleSolution_adjoint
end

ZygoteRules.@adjoint function DiffEqBase.EnsembleSolution(sim, time, converged)
    out = EnsembleSolution(sim, time, converged)
    function EnsembleSolution_adjoint(p̄::AbstractArray{T, N}) where {T, N}
        arrarr = [[p̄[ntuple(x -> Colon(), Val(N - 2))..., j, i]
                   for j in 1:size(p̄)[end - 1]] for i in 1:size(p̄)[end]]
        (EnsembleSolution(arrarr, 0.0, true), nothing, nothing)
    end
    function EnsembleSolution_adjoint(p̄::AbstractArray{<:AbstractArray, 1})
        (EnsembleSolution(p̄, 0.0, true), nothing, nothing)
    end
    function EnsembleSolution_adjoint(p̄::EnsembleSolution)
        (p̄, nothing, nothing)
    end
    out, EnsembleSolution_adjoint
end

ZygoteRules.@adjoint function ZygoteRules.literal_getproperty(sim::EnsembleSolution,
                                                              ::Val{:u})
    sim.u, p̄ -> (EnsembleSolution(p̄, 0.0, true),)
end

#=
ChainRulesCore.frule(f::ODEFunction,u,p,t)
  if f.jvp === nothing
    ChainRulesCore.frule(f.f,u,p,t)
  else
    function ode_jvp(f,du,dp,dt)
      f.jvp_u(du,u,p,t) + f.jvp_p(dp,u,p,t) + f.tgrad(u,p,t)*dt
    end
    f.f(u,p,t),ode_jvp
end
=#

#=
function ChainRulesCore.rrule(f::ODEFunction,u,p,t)
  if f.vjp === nothing
    ChainRulesCore.rrule(f.f,u,p,t)
  else
    f.vjp(u,p,t)
  end
end
=#

#=
ZygoteRules.@adjoint function (f::ODEFunction)(u,p,t)
  if f.vjp === nothing
    ZygoteRules._pullback(f.f,u,p,t)
  else
    f.vjp(u,p,t)
  end
end
=#

#=
ZygoteRules.@adjoint! function (f::ODEFunction)(du,u,p,t)
  if f.vjp === nothing
    ZygoteRules._pullback(f.f,du,u,p,t)
  else
    f.vjp(du,u,p,t)
  end
end
=#

ChainRulesCore.rrule(::typeof(numargs), f) = (numargs(f), df -> (NoTangent(), NoTangent()))
