#=
ZygoteRules.@adjoint function ODESolution(u,args...)
  function ODESolutionAdjoint(ȳ)
    (ȳ,ntuple(_->nothing, length(args))...)
  end
  ODESolution(u,args...), ODESolutionAdjoint
end
=#

ZygoteRules.@adjoint function ODESolution{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11
                }(u,args...) where {T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11}

                function ODESolutionAdjoint(ȳ)
                  (ȳ,ntuple(_->nothing, length(args))...)
                end

                ODESolution{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11}(u,args...),ODESolutionAdjoint
end

ZygoteRules.@adjoint function ZygoteRules.literal_getproperty(sol::AbstractTimeseriesSolution, ::Val{:u})
  function solu_adjoint(Δ)
        zerou = zero(sol.prob.u0)
        _Δ = @. ifelse(Δ == nothing,(zerou,),Δ)
        (DiffEqBase.build_solution(sol.prob,sol.alg,sol.t,_Δ),)
  end
  sol.u,solu_adjoint
end

ZygoteRules.@adjoint function ZygoteRules.literal_getproperty(sol::AbstractNoTimeSolution, ::Val{:u})
  function solu_adjoint(Δ)
        zerou = zero(sol.prob.u0)
        _Δ = @. ifelse(Δ == nothing,(zerou,),Δ)
        (DiffEqBase.build_solution(sol.prob,sol.alg,_Δ,sol.resid),)
  end
  sol.u,solu_adjoint
end

ZygoteRules.@adjoint function (f::ODEFunction)(u,p,t)
  if f.vjp === nothing
    ZygoteRules._pullback(f.f,u,p,t)
  else
    f.vjp(u,p,t)
  end
end

ZygoteRules.@adjoint! function (f::ODEFunction)(du,u,p,t)
  if f.vjp === nothing
    ZygoteRules._pullback(f.f,du,u,p,t)
  else
    f.vjp(du,u,p,t)
  end
end

ZygoteRules.@adjoint function SciMLBase.tmap(f, args::Union{AbstractArray,Tuple}...)
  ∇tmap(__context__, f, args...)
end

ZygoteRules.@adjoint function DiffEqBase.EnsembleSolution(sim,time,converged)
  out = EnsembleSolution(sim,time,converged)
  function EnsembleSolution_adjoint(p̄::AbstractArray{T,N}) where {T,N}
    arrarr = [[p̄[ntuple(x->Colon(),Val(N-2))...,j,i] for j in 1:size(p̄)[end-1]] for i in 1:size(p̄)[end]]
    (EnsembleSolution(arrarr, 0.0, true),nothing,nothing)
  end
  function EnsembleSolution_adjoint(p̄::EnsembleSolution)
    (p̄,nothing,nothing)
  end
  out,EnsembleSolution_adjoint
end
ZygoteRules.@adjoint ZygoteRules.literal_getproperty(sim::EnsembleSolution, ::Val{:u}) = sim.u, p̄ -> (EnsembleSolution(p̄, 0.0, true),)

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

function ChainRulesCore.rrule(f::ODEFunction,u,p,t)
  if f.vjp === nothing
    ChainRulesCore.rrule(f.f,u,p,t)
  else
    f.vjp(u,p,t)
  end
end

ZygoteRules.@adjoint numargs(f) = (numargs(f),df->(nothing,))
ChainRulesCore.rrule(::typeof(numargs),f) = (numargs(f),df->(nothing,))

# Until https://github.com/FluxML/Zygote.jl/issues/664 is fixed
ZygoteRules.@adjoint function Base.pairs(x::T) where T
         y = Base.pairs(x)
         back(dx::NamedTuple) = (dx.data,)
         function back(dx::Dict)
           T <: AbstractDict && return (dx,)
           z = zero(x)
           for (k,v) in dx
             z[k] = v
           end
           (z,)
         end
         y, back
       end
