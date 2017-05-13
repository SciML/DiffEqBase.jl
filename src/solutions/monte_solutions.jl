type MonteCarloTestSolution{T,N,S} <: AbstractMonteCarloSolution{T,N}
  u::S
  errors::Dict{Symbol,Vector{T}}
  error_means::Dict{Symbol,T}
  error_medians::Dict{Symbol,T}
  elapsedTime::Float64
  converged::Bool
end
function MonteCarloTestSolution{T,N}(sim::AbstractMonteCarloSolution{T,N},errors,error_means,error_medians)
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

MonteCarloSolution{N}(sim, dims::NTuple{N},elapsedTime) =
                 MonteCarloSolution{eltype(eltype(sim)), N, typeof(sim)}(sim,elapsedTime,false)
MonteCarloSolution(sim,elapsedTime,converged) =
            MonteCarloSolution(sim, (size(sim[1])..., length(sim)),elapsedTime,false)

type MonteCarloSummary{T,N,Tt,S,S2} <: AbstractMonteCarloSolution{T,N}
  t::Tt
  u::S
  v::S2
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

@recipe function f(sim::AbstractMonteCarloSolution;idxs=eachindex(sim.u[1]))
  for i in idxs
    @series begin
      legend := false
      sol[i]
    end
  end
end

@recipe function f(sim::MonteCarloSummary;
                   idxs= typeof(sim.u[1])<:AbstractArray ? eachindex(sim.u[1]) : 1,
                   error_style=:ribbon)
  if typeof(sim.u[1])<:AbstractArray
    ci = vecarr_to_vectors(VectorOfArray([sqrt(sim.v[i]).*1.96 for i in 1:length(sim.v)]))
  else
    ci = [[sqrt(sim.v[i]).*1.96 for i in 1:length(sim.v)]]
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
        ribbon --> ci[i]
      elseif error_style == :bars
        yerr --> ci[i]
      elseif error_style == :none
        nothing
      else
        error("error_style not recognized")
      end
      sim.t,u[i]
    end
  end
end
