struct MonteCarloTestSolution{T,N,S} <: AbstractMonteCarloSolution{T,N}
  u::S
  errors::Dict{Symbol,Vector{T}}
  weak_errors::Dict{Symbol,T}
  error_means::Dict{Symbol,T}
  error_medians::Dict{Symbol,T}
  elapsedTime::Float64
  converged::Bool
end
function MonteCarloTestSolution(sim::AbstractMonteCarloSolution{T,N},errors,weak_errors,error_means,error_medians,elapsedTime,converged) where {T,N}
  MonteCarloTestSolution{T,N,typeof(sim.u)}(sim.u,errors,weak_errors,error_means,error_medians,sim.elapsedTime,sim.converged)
end
function MonteCarloTestSolution(u,errors,weak_errors,error_means,error_medians,elapsedTime,converged)
  MonteCarloTestSolution(MonteCarloSolution(u,elapsedTime,converged),errors,weak_errors,error_means,error_medians,elapsedTime,converged)
end

struct MonteCarloSolution{T,N,S} <: AbstractMonteCarloSolution{T,N}
  u::S
  elapsedTime::Float64
  converged::Bool
end
MonteCarloSolution(sim, dims::NTuple{N},elapsedTime,converged) where {N} =
                  MonteCarloSolution{eltype(eltype(sim)), N, typeof(sim)}(sim,elapsedTime,converged)
MonteCarloSolution(sim,elapsedTime,converged) =
             MonteCarloSolution(sim, (length(sim),),elapsedTime,converged) # Vector of some type which is not an array
MonteCarloSolution(sim::T,elapsedTime,converged) where T <: AbstractVector{T2} where T2 <: AbstractArray =
             MonteCarloSolution(sim, (size(sim[1])..., length(sim)),elapsedTime,converged) # Requires `size` defined on `sim`

struct MonteCarloSummary{T,N,Tt,S,S2,S3,S4} <: AbstractMonteCarloSolution{T,N}
  t::Tt
  u::S
  v::S2
  qlow::S3
  qhigh::S4
  num_monte::Int
  elapsedTime::Float64
  converged::Bool
end

function calculate_monte_errors(sim::AbstractMonteCarloSolution;kwargs...)
  calculate_monte_errors(sim.u;elapsedTime=sim.elapsedTime,converged=sim.converged,kwargs...)
end

function calculate_monte_errors(u;elapsedTime=0.0,converged=false,
                                weak_timeseries_errors=false,weak_dense_errors=false)
  errors = Dict{Symbol,Vector{eltype(u[1].u[1])}}() #Should add type information
  error_means  = Dict{Symbol,eltype(u[1].u[1])}()
  error_medians= Dict{Symbol,eltype(u[1].u[1])}()
  for k in keys(u[1].errors)
    errors[k] = [sol.errors[k] for sol in u]
    error_means[k] = mean(errors[k])
    error_medians[k]=median(errors[k])
  end
  # Now Calculate Weak Errors
  weak_errors =  Dict{Symbol,eltype(u[1].u[1])}()
  # Final
  m_final = mean([s[end] for s in u])
  m_final_analytic = mean([s.u_analytic[end] for s in u])
  res = norm(m_final - m_final_analytic)
  weak_errors[:weak_final] = res
  if weak_timeseries_errors
    ts_weak_errors = [mean([u[j][i] - u[j].u_analytic[i] for j in 1:length(u)]) for i in 1:length(u[1])]
    ts_l2_errors = [sqrt.(sum(abs2,err)/length(err)) for err in ts_weak_errors]
    l2_tmp = sqrt(sum(abs2,ts_l2_errors)/length(ts_l2_errors))
    max_tmp = maximum([maximum(abs.(err)) for err in ts_weak_errors])
    weak_errors[:weak_l2] = l2_tmp
    weak_errors[:weak_l∞] = max_tmp
  end
  if weak_dense_errors
    densetimes = collect(range(u[1].t[1], stop=u[1].t[end], length=100))
    u_analytic = [[sol.prob.f(Val{:analytic},sol.prob.u0,sol.prob.p,densetimes[i],sol.W(densetimes[i])[1]) for i in eachindex(densetimes)] for sol in u]
    dense_weak_errors = [mean([u[j](densetimes)[i] - u_analytic[j][i] for j in 1:length(u)]) for i in eachindex(densetimes)]
    dense_L2_errors = [sqrt.(sum(abs2,err)/length(err)) for err in dense_weak_errors]
    L2_tmp = sqrt(sum(abs2,dense_L2_errors)/length(dense_L2_errors))
    max_tmp = maximum([maximum(abs.(err)) for err in dense_weak_errors])
    weak_errors[:weak_L2] = L2_tmp
    weak_errors[:weak_L∞] = max_tmp
  end
  return MonteCarloTestSolution(u,errors,weak_errors,error_means,error_medians,elapsedTime,converged)
end

### Displays

Base.summary(A::AbstractMonteCarloSolution) = string("MonteCarloSolution Solution of length ",length(A.u)," with uType:\n",eltype(A.u))
function Base.show(io::IO, A::AbstractMonteCarloSolution)
  print(io,summary(A))
end
function Base.show(io::IO, m::MIME"text/plain", A::AbstractMonteCarloSolution)
  print(io,summary(A))
end

### Plot Recipes

@recipe function f(sim::AbstractMonteCarloSolution;
                   zcolors = typeof(sim.u)<:AbstractArray ? fill(nothing, length(sim.u)) : nothing,
                   idxs = typeof(sim.u)<:AbstractArray ? eachindex(sim.u) : 1)
  for i in idxs
    size(sim[i].u, 1) == 0 && continue
    @series begin
      legend := false
      xlims --> (-Inf,Inf)
      ylims --> (-Inf,Inf)
      zlims --> (-Inf,Inf)
      zcolor --> zcolors[i]
      sim[i]
    end
  end
end

@recipe function f(sim::MonteCarloSummary;
                   idxs= typeof(sim.u[1])<:AbstractArray ? eachindex(sim.u[1]) : 1,
                   error_style=:ribbon,ci_type=:quantile)
  if typeof(sim.u[1])<:AbstractArray
    u = vecarr_to_vectors(sim.u)
  else
    u = [sim.u.u]
  end
  if ci_type == :SEM
    if typeof(sim.u[1])<:AbstractArray
      ci_low = vecarr_to_vectors(VectorOfArray([sqrt(sim.v[i]/sim.num_monte).*1.96 for i in 1:length(sim.v)]))
      ci_high = ci_low
    else
      ci_low = [[sqrt(sim.v[i]/length(sim.num_monte)).*1.96 for i in 1:length(sim.v)]]
      ci_high = ci_low
    end
  elseif ci_type == :quantile
    if typeof(sim.u[1])<:AbstractArray
      ci_low = u - vecarr_to_vectors(sim.qlow)
      ci_high = vecarr_to_vectors(sim.qhigh) - u
    else
      ci_low = [u[1] - sim.qlow.u]
      ci_high = [sim.qhigh.u - u[1]]
    end
  else
    error("ci_type choice not valid. Must be :variance or :quantile")
  end
  for i in idxs
    @series begin
      legend --> false
      lw --> 3
      fillalpha --> 0.2
      if error_style == :ribbon
        ribbon --> (ci_low[i],ci_high[i])
      elseif error_style == :bars
        yerror --> (ci_low[i],ci_high[i])
      elseif error_style == :none
        nothing
      else
        error("error_style not recognized")
      end
      sim.t,u[i]
    end
  end
end
