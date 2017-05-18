type MonteCarloTestSolution{T,N,S} <: AbstractMonteCarloSolution{T,N}
  u::S
  errors::Dict{Symbol,Vector{T}}
  error_means::Dict{Symbol,T}
  error_medians::Dict{Symbol,T}
  elapsedTime::Float64
  converged::Bool
end
function MonteCarloTestSolution{T,N}(sim::AbstractMonteCarloSolution{T,N},errors,error_means,error_medians,converged)
  MonteCarloTestSolution{T,N,typeof(sim.u)}(sim.u,errors,error_means,error_medians,sim.elapsedTime,sim.converged)
end

type MonteCarloSolution{T,N,S} <: AbstractMonteCarloSolution{T,N}
  u::S
  elapsedTime::Float64
  converged::Bool
end
MonteCarloSolution{N}(sim, dims::NTuple{N},elapsedTime,converged) =
                  MonteCarloSolution{eltype(eltype(sim)), N, typeof(sim)}(sim,elapsedTime,converged)
MonteCarloSolution(sim,elapsedTime,converged) =
             MonteCarloSolution(sim, (size(sim[1])..., length(sim)),elapsedTime,converged)

type MonteCarloSummary{T,N,Tt,S,S2,S3,S4} <: AbstractMonteCarloSolution{T,N}
  t::Tt
  u::S
  v::S2
  qlow::S3
  qhigh::S4
  num_monte::Int
  elapsedTime::Float64
  converged::Bool
end

function calculate_monte_errors(sim::AbstractMonteCarloSolution)
  u = sim.u
  errors = Dict{Symbol,Vector{eltype(u[1].u[1])}}() #Should add type information
  error_means  = Dict{Symbol,eltype(u[1].u[1])}()
  error_medians= Dict{Symbol,eltype(u[1].u[1])}()
  for k in keys(u[1].errors)
    errors[k] = [sol.errors[k] for sol in u]
    error_means[k] = mean(errors[k])
    error_medians[k]=median(errors[k])
  end
  return MonteCarloTestSolution(sim,errors,error_means,error_medians,sim.converged)
end

@recipe function f(sim::AbstractMonteCarloSolution;
                   idxs = typeof(sim.u)<:AbstractArray ? eachindex(sim.u) : 1)
  for i in idxs
    @series begin
      legend := false
      sim[i]
    end
  end
end

@recipe function f(sim::MonteCarloSummary;
                   idxs= typeof(sim.u[1])<:AbstractArray ? eachindex(sim.u[1]) : 1,
                   error_style=:ribbon,ci_type=:quantile)
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
      ci_low = vecarr_to_vectors(sim.qlow)
      ci_high = vecarr_to_vectors(sim.qhigh)
    else
      ci_low = [sim.qlow]
      ci_high = [sim.qhigh]
    end
  else
    error("ci_type choice not valid. Must be :variance or :quantile")
  end
  if typeof(sim.u[1])<:AbstractArray
    u = vecarr_to_vectors(sim.u)
  else
    u = [sim.u.u]
  end
  for i in idxs
    @series begin
      legend := false
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
