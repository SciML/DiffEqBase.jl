type MonteCarloTestSolution{S,T} <: AbstractMonteCarloSolution
  solution_data::S
  errors::Dict{Symbol,Vector{T}}
  error_means::Dict{Symbol,T}
  error_medians::Dict{Symbol,T}
  elapsedTime::Float64
end

type MonteCarloSolution{T} <: AbstractMonteCarloSolution
  solution_data::Vector{T}
  elapsedTime::Float64
end

function calculate_monte_errors(sim::AbstractMonteCarloSolution)
  solution_data = sim.solution_data
  errors = Dict{Symbol,Vector{eltype(solution_data[1].u[1])}}() #Should add type information
  error_means  = Dict{Symbol,eltype(solution_data[1].u[1])}()
  error_medians= Dict{Symbol,eltype(solution_data[1].u[1])}()
  for k in keys(solution_data[1].errors)
    errors[k] = [sol.errors[k] for sol in solution_data]
    error_means[k] = mean(errors[k])
    error_medians[k]=median(errors[k])
  end
  return MonteCarloTestSolution(solution_data,errors,error_means,error_medians,sim.elapsedTime)
end

Base.length(sim::AbstractMonteCarloSolution) = length(sim.solution_data)
Base.endof( sim::AbstractMonteCarloSolution) = length(sim)
Base.getindex(sim::AbstractMonteCarloSolution,i::Int) = sim.solution_data[i]
Base.getindex(sim::AbstractMonteCarloSolution,i::Int,I::Int...) = sim.solution_data[i][I...]
Base.size(sim::AbstractMonteCarloSolution) = (length(sim),)
Base.start(sim::AbstractMonteCarloSolution) = 1
function Base.next(sim::AbstractMonteCarloSolution,state)
  state += 1
  (sim[state],state)
end
Base.done(sim::AbstractMonteCarloSolution,state) = state >= length(sim)

@recipe function f(sim::AbstractMonteCarloSolution)

  for sol in sim
    @series begin
      legend := false
      sol
    end
  end
end
