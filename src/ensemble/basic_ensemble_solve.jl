"""
$(TYPEDEF)
"""
abstract type BasicEnsembleAlgorithm <: EnsembleAlgorithm end

"""
$(TYPEDEF)
"""
struct EnsembleThreads <: BasicEnsembleAlgorithm end

"""
$(TYPEDEF)
"""
struct EnsembleDistributed <: BasicEnsembleAlgorithm end

"""
$(TYPEDEF)
"""
struct EnsembleSplitThreads <: BasicEnsembleAlgorithm end

"""
$(TYPEDEF)
"""
struct EnsembleSerial <: BasicEnsembleAlgorithm end

#=
if (kwargs[:parallel_type] == != :none && kwargs[:parallel_type] == != :threads)
  error("Distributed arrays cannot be generated via none or threads")
end
(batch_size != trajectories) && warn("batch_size and reductions are ignored when !collect_result")

elapsed_time = @elapsed u = DArray((trajectories,)) do II
    solve_batch(prob,alg,kwargs[:parallel_type] ==,II[1],pmap_batch_size,kwargs...)
end
return EnsembleSolution(u,elapsed_time,false)
=#

function __solve(prob::AbstractEnsembleProblem,
                 alg::Union{DEAlgorithm,Nothing};
                 kwargs...)
    #=
    if alg isa EnsembleAlgorithm
      @error "You forgot to pass a DE solver algorithm! Only a EnsembleAlgorithm has been supplied. Exiting"
    end
    =#
    if :parallel_type ∈ keys(kwargs)
      if kwargs[:parallel_type] == :none
        ensemblealg = EnsembleSerial()
      elseif kwargs[:parallel_type] == :pmap || kwargs[:parallel_type] == :parfor
        ensemblealg = EnsembleDistributed()
      elseif kwargs[:parallel_type] == :threads
        ensemblealg = EnsembleThreads()
      elseif kwargs[:parallel_type] == :split_threads
        ensemblealg = EnsembleSplitThreads()
      else
        @error "parallel_type value not recognized"
      end
    elseif alg isa EnsembleAlgorithm
      # Assume DifferentialEquations.jl is being used, so default alg
      ensemblealg = alg
      alg = nothing
    else
      ensemblealg = EnsembleThreads()
    end
    if :num_monte ∈ keys(kwargs)
      @warn "num_monte has been replaced by trajectories"
      trajectories = kwargs[:num_monte]
    else
      @assert :trajectories ∈ keys(kwargs)
      trajectories = kwargs[:trajectories]
    end
    if :parallel_type ∈ keys(kwargs)
      @warn "parallel_type has been deprecated. Please refer to the docs for the new dispatch-based system."
    end
    __solve(prob,alg,ensemblealg;trajectories=trajectories,kwargs...)
end

tighten_container_eltype(u::Vector{Any}) = map(identity, u)
tighten_container_eltype(u) = u

function __solve(prob::AbstractEnsembleProblem,
                 alg::Union{DEAlgorithm,Nothing},
                 ensemblealg::BasicEnsembleAlgorithm;
                 trajectories, batch_size = trajectories,
                 pmap_batch_size = batch_size÷100 > 0 ? batch_size÷100 : 1, kwargs...)

  num_batches = trajectories ÷ batch_size
  num_batches < 1 && error("trajectories ÷ batch_size cannot be less than 1, got $num_batches")
  num_batches * batch_size != trajectories && (num_batches += 1)

  function batch_function(II)
    batch_data = solve_batch(prob,alg,ensemblealg,II,pmap_batch_size;kwargs...)
  end

  if num_batches == 1 && prob.reduction === DEFAULT_REDUCTION
    elapsed_time = @elapsed u = batch_function(1:trajectories)
    _u = tighten_container_eltype(u)
    return EnsembleSolution(_u,elapsed_time,true)
  end

  converged::Bool = false
  i = 1
  II = (batch_size*(i-1)+1):batch_size*i

  batch_data = batch_function(II)
  u = prob.u_init === nothing ? similar(batch_data, 0) : prob.u_init
  u,converged = prob.reduction(u,batch_data,II)
  elapsed_time = @elapsed for i in 2:num_batches
    converged && break
    if i == num_batches
      II = (batch_size*(i-1)+1):trajectories
    else
      II = (batch_size*(i-1)+1):batch_size*i
    end
    batch_data = batch_function(II)
    u,converged = prob.reduction(u,batch_data,II)
  end

  u = reduce(vcat, u)
  _u = tighten_container_eltype(u)

  return EnsembleSolution(_u,elapsed_time,converged)

end

function batch_func(i,prob,alg;kwargs...)
  iter = 1
  _prob = prob.safetycopy ? deepcopy(prob.prob) : prob.prob
  new_prob = prob.prob_func(_prob,i,iter)
  rerun = true
  x = prob.output_func(solve(new_prob,alg;kwargs...),i)
  if !(typeof(x) <: Tuple)
      @warn("output_func should return (out,rerun). See docs for updated details")
      _x = (x,false)
  else
    _x = x
  end
  rerun = _x[2]
  while rerun
      iter += 1
      _prob = prob.safetycopy ? deepcopy(prob.prob) : prob.prob
      new_prob = prob.prob_func(_prob,i,iter)
      x = prob.output_func(solve(new_prob,alg;kwargs...),i)
      if !(typeof(x) <: Tuple)
          @warn("output_func should return (out,rerun). See docs for updated details")
          _x = (x,false)
      else
        _x = x
      end
      rerun = _x[2]
  end
  _x[1]
end

function solve_batch(prob,alg,ensemblealg::EnsembleDistributed,II,pmap_batch_size;kwargs...)
  wp=CachingPool(workers())
  batch_data = pmap(wp,II,batch_size=pmap_batch_size) do i
    batch_func(i,prob,alg;kwargs...)
  end
  map(identity,batch_data)
end

function solve_batch(prob,alg,::EnsembleSerial,II,pmap_batch_size;kwargs...)
  batch_data = map(II) do i
    batch_func(i,prob,alg;kwargs...)
  end
  batch_data
end

function solve_batch(prob,alg,ensemblealg::EnsembleThreads,II,pmap_batch_size;kwargs...)
  function multithreaded_batch(batch_idx)
    i = II[batch_idx]
    iter = 1
    _prob = prob.safetycopy ? deepcopy(prob.prob) : prob.prob
    new_prob = prob.prob_func(_prob,i,iter)
    x = prob.output_func(solve(new_prob,alg;kwargs...),i)
    if !(typeof(x) <: Tuple)
        @warn("output_func should return (out,rerun). See docs for updated details")
        _x = (x,false)
    else
      _x = x
    end
    rerun = _x[2]

    while rerun
        iter += 1
        _prob = prob.safetycopy ? deepcopy(prob.prob) : prob.prob
        new_prob = prob.prob_func(_prob,i,iter)
        x = prob.output_func(solve(new_prob,alg;kwargs...),i)
        if !(typeof(x) <: Tuple)
            @warn("output_func should return (out,rerun). See docs for updated details")
            _x = (x,false)
        else
          _x = x
        end
        rerun = _x[2]
    end
    _x[1]
  end

  batch_data = Vector{Core.Compiler.return_type(multithreaded_batch,Tuple{typeof(first(II))})}(undef,length(II))
  let
    Threads.@threads for batch_idx in axes(batch_data, 1)
        batch_data[batch_idx] = multithreaded_batch(batch_idx)
    end
  end
  batch_data
end

function solve_batch(prob,alg,::EnsembleSplitThreads,II,pmap_batch_size;kwargs...)
  wp=CachingPool(workers())
  N = nworkers()
  batch_size = length(II)÷N
  batch_data = let
    pmap(wp,1:N,batch_size=pmap_batch_size) do i
      if i == N
        I_local = II[(batch_size*(i-1)+1):end]
      else
        I_local = II[(batch_size*(i-1)+1):(batch_size*i)]
      end
      thread_monte(prob,I_local,alg,i;kwargs...)
    end
  end
  reduce(vcat,batch_data)
end

function thread_monte(prob,II,alg,procid;kwargs...)
  function multithreaded_batch(j)
    i = II[j]
    iter = 1
    _prob = prob.safetycopy ? deepcopy(prob.prob) : prob.prob
    new_prob = prob.prob_func(_prob,i,iter)
    rerun = true
    x = prob.output_func(solve(new_prob,alg;kwargs...),i)
    if !(typeof(x) <: Tuple)
        @warn("output_func should return (out,rerun). See docs for updated details")
        _x = (x,false)
    else
      _x = x
    end
    rerun = _x[2]
    while rerun
        iter += 1
        _prob = prob.safetycopy ? deepcopy(prob.prob) : prob.prob
        new_prob = prob.prob_func(_prob,i,iter)
        x = prob.output_func(solve(new_prob,alg;kwargs...),i)
        if !(typeof(x) <: Tuple)
            @warn("output_func should return (out,rerun). See docs for updated details")
            _x = (x,false)
        else
          _x = x
        end
        rerun = _x[2]
    end
    _x[1]
  end

  batch_data = Vector{Core.Compiler.return_type(multithreaded_batch,Tuple{typeof(first(II))})}(undef,length(II))

  let
    Threads.@threads for j in 1:length(II)
      batch_data[j] = multithreaded_batch(j)
    end
  end
  batch_data
end
