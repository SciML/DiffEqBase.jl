## preamble!

@muladd function initialize_cache!(nlcache::Union{NLFunctionalConstantCache,NLFunctionalCache},
                                   nlsolver::NLSolver{<:NLFunctional}, integrator)
  nlcache.tstep = integrator.t + nlsolver.c * integrator.dt

  nothing
end

@muladd function initialize_cache!(nlcache::Union{NLAndersonConstantCache,NLAndersonCache},
                                   nlsolver::NLSolver{<:NLAnderson}, integrator)
  nlcache.history = 0
  nlcache.tstep = integrator.t + nlsolver.c * integrator.dt

  nothing
end

## apply_step!

function apply_step!(nlsolver::NLSolver{<:NLAnderson,false}, integrator)
  # perform Anderson acceleration
  nlsolver.gz = anderson(nlsolver.gz, nlsolver.iter, nlsolver.alg.aa_start, nlsolver.cache, integrator)

  # apply step
  _apply_step!(nlsolver, integrator)

  nothing
end

function apply_step!(nlsolver::NLSolver{<:NLAnderson,true}, integrator)
  # perform Anderson acceleration
  anderson!(nlsolver.gz, nlsolver.iter, nlsolver.alg.aa_start, nlsolver.cache, integrator)

  # apply step
  _apply_step!(nlsolver, integrator)

  nothing  
end

## perform_step

"""
    perform_step!(nlsolver::NLSolver{<:Union{NLFunctional,NLAnderson}}, integrator)

Update `nlsolver.z` with the next iterate `g(nlsolver.z)`, where
```math
g(z) = dt⋅f(tmp + γ⋅z, p, t + c⋅dt).
```
"""
@muladd function perform_step!(nlsolver::NLSolver{<:Union{NLFunctional,NLAnderson},false}, integrator)
  @unpack p,dt = integrator
  @unpack z,tmp,γ,cache = nlsolver
  @unpack tstep = cache

  mass_matrix = integrator.f.mass_matrix
  f = nlsolve_f(integrator)

  # evaluate function
  ztmp = @.. tmp + γ * z
  if mass_matrix == I
    gz = dt .* f(ztmp, p, tstep)
    dz = gz .- z
  else
    mz = _reshape(mass_matrix * _vec(z), axes(z))
    dz = dt .* f(ztmp, p, tstep) .- mz
    gz = z .+ dz
  end
  if has_destats(integrator)
    integrator.destats.nf += 1
  end

  # cache residuals 
  cache.dz = dz

  # save next step
  nlsolver.gz = gz

  nothing
end

@muladd function perform_step!(nlsolver::NLSolver{<:Union{NLFunctional,NLAnderson},true}, integrator)
  @unpack p,dt = integrator
  @unpack z,gz,tmp,γ,cache = nlsolver
  @unpack dz,tstep,k = cache

  mass_matrix = integrator.f.mass_matrix
  f = nlsolve_f(integrator)

  # use gz as temporary variable
  ztmp = gz

  # evaluate function
  @.. ztmp = tmp + γ * z
  f(k, ztmp, p, tstep)
  if has_destats(integrator)
    integrator.destats.nf += 1
  end
  
  if mass_matrix == I
    @.. gz = dt * k
    @.. dz = gz - z
  else
    mul!(vec(ztmp), mass_matrix, vec(z))
    @.. dz = dt * k - ztmp
    @.. gz = z + dz
  end

  nothing
end

## resize!

function Base.resize!(nlcache::NLFunctionalCache, i::Int)
  resize!(nlcache.dz, i)
  resize!(nlcache.k, i)
  resize!(nlcache.atmp, i)
  nothing
end

function Base.resize!(nlcache::Union{NLAndersonCache,NLAndersonConstantCache}, nlsolver::NLSolver{<:NLAnderson},
  integrator, i::Int)
  resize!(nlcache, nlsolver.alg, i)
end

function Base.resize!(nlcache::NLAndersonCache, nlalg::NLAnderson, i::Int)
  resize!(nlcache.dz, i)
  resize!(nlcache.k, i)
  resize!(nlcache.atmp, i)
  resize!(nlcache.gzprev, i)
  resize!(nlcache.dzprev, i)

  # determine new maximum history
  max_history_old = length(Δgzs)
  max_history = min(nlalg.max_history, nlalg.maxiters, i)

  resize!(nlcache.γs, max_history)
  resize!(nlcache.Δgzs, max_history)
  if max_history > max_history_old
    for i in (max_history_old + 1):max_history
      nlcache.Δgzs[i] = zero(nlcache.gzprev)
    end
  end

  if max_history != max_history_old
    nlcache.Q = typeof(nlcache.Q)(undef, i, max_history)
    nlcache.R = typeof(nlcache.R)(undef, max_history, max_history)
  end

  nothing
end

function Base.resize!(nlcache::NLAndersonConstantCache, nlalg::NLAnderson, i::Int)
  # determine new maximum history
  max_history_old = length(nlcache.Δgzs)
  max_history = min(nlalg.max_history, nlalg.maxiters, i)

  resize!(nlcache.γs, max_history)
  resize!(nlcache.Δgzs, max_history)

  if max_history != max_history_old
    nlcache.Q = typeof(nlcache.Q)(undef, i, max_history)
    nlcache.R = typeof(nlcache.R)(undef, max_history, max_history)
  end

  nothing
end