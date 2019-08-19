## preamble!

@muladd function preamble!(nlsolver::NLSolver{<:NLFunctional}, integrator)
  nlsolver.cache.tstep = integrator.t + nlsolver.c * integrator.dt

  nothing
end

@muladd function preamble!(nlsolver::NLSolver{<:NLAnderson}, integrator)
  @unpack cache = nlsolver

  cache.history = 0
  cache.tstep = integrator.t + nlsolver.c * integrator.dt

  nothing
end

## loopheader!

@muladd function loopheader!(nlsolver::NLSolver{<:NLAnderson,false}, integrator, iter::Int)
  @unpack aa_start = nlsolver.alg
  
  # perform Anderson acceleration
  if iter == aa_start
    @unpack cache = nlsolver

    # update cached values for next step of Anderson acceleration
    cache.dzprev = cache.dz
    cache.gprev = nlsolver.z
  elseif iter > aa_start
    @unpack z,cache = nlsolver
    @unpack dz,Δgs,Q,R,γs,history,droptol = cache

    # increase size of history
    history += 1
  
    # remove oldest history if maximum size is exceeded
    max_history = length(Δgs)
    if history > max_history
      # circularly shift differences of G(z)
      for i in 1:(max_history-1)
        Δgs[i] = Δgs[i + 1]
      end
  
      # delete left-most column of QR decomposition
      qrdelete!(Q, R, max_history)
  
      # update size of history
      history = max_history
    end
  
    # update history of differences of G(z)
    Δgs[history] = @.. z - cache.gprev
  
    # replace/add difference of residuals as right-most column to QR decomposition
    qradd!(Q, R, _vec(dz .- cache.dzprev), history)
  
    # update cached values
    cache.dzprev = dz
    cache.gprev = z
  
    # define current Q and R matrices
    Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))
  
    # check condition (TODO: incremental estimation)
    if droptol !== nothing
      while cond(R) > droptol && history > 1
        qrdelete!(Q, R, history)
        history -= 1
        Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))
      end
    end
  
    # solve least squares problem
    γscur = view(γs, 1:history)
    ldiv!(Rcur, mul!(γscur, Qcur', _vec(dz)))
    if has_destats(integrator)
      integrator.destats.nsolve += 1
    end
  
    # update iterate
    for i in 1:history
      z = @.. z - γs[i] * Δgs[i]
    end
  
    # update cached values
    nlsolver.z = z
    cache.history = history
  end

  # apply step
  apply_step!(nlsolver, integrator)
end

@muladd function loopheader!(nlsolver::NLSolver{<:NLAnderson,true}, integrator, iter::Int)
  @unpack aa_start = nlsolver.alg

  # perform Anderson acceleration
  if iter == aa_start
    @unpack z,cache = nlsolver

    # update cached values for next step of Anderson acceleration
    @.. cache.gprev = z
    @.. cache.dzprev = cache.dz
  elseif iter > aa_start
    @unpack z,cache = nlsolver
    @unpack gprev,dz,dzprev,Δgs,Q,R,γs,history,droptol = cache

    # increase size of history
    history += 1

    # remove oldest history if maximum size is exceeded
    max_history = length(Δgs)
    if history > max_history
      # circularly shift differences of z
      ptr = Δgs[1]
      for i in 1:(max_history-1)
        Δgs[i] = Δgs[i + 1]
      end
      Δgs[max_history] = ptr

      # delete left-most column of QR decomposition
      qrdelete!(Q, R, max_history)

      # update size of history
      history = max_history
    end

    # update history of differences of z
    @.. Δgs[history] = z - gprev

    # replace/add difference of residuals as right-most column to QR decomposition
    @.. dzprev = dz - dzprev
    qradd!(Q, R, vec(dzprev), history)

    # update cached values
    @.. dzprev = dz
    @.. gprev = z

    # define current Q and R matrices
    Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))

    # check condition (TODO: incremental estimation)
    if droptol !== nothing
      while cond(R) > droptol && history > 1
        qrdelete!(Q, R, history)
        history -= 1
        Qcur, Rcur = view(Q, :, 1:history), UpperTriangular(view(R, 1:history, 1:history))
      end
    end

    # solve least squares problem
    γscur = view(γs, 1:history)
    ldiv!(Rcur, mul!(γscur, Qcur', vec(dz)))
    if has_destats(integrator)
      integrator.destats.nsolve += 1
    end

    # update next iterate
    for i in 1:history
      @.. z = z - γs[i] * Δgs[i]
    end

    # update cached values
    cache.history = history
  end

  # apply step
  apply_step!(nlsolver, integrator)
end

## perform_step

"""
    perform_step!(nlsolver::NLSolver{<:Union{NLFunctional,NLAnderson}}, integrator, iter::Int)

Update `nlsolver.z` with the next iterate `G(nlsolver.z)`, where
```math
G(z) = dt⋅f(tmp + γ⋅z, p, t + c⋅dt).
```
"""
@muladd function perform_step!(nlsolver::NLSolver{<:Union{NLFunctional,NLAnderson},false}, integrator, iter::Int)
  @unpack p,dt = integrator
  @unpack zprev,tmp,γ,cache = nlsolver
  @unpack tstep = cache

  mass_matrix = integrator.f.mass_matrix
  f = nlsolve_f(integrator)

  # evaluate function
  ztmp = @.. tmp + γ * zprev
  if mass_matrix == I
    z = dt .* f(ztmp, p, tstep)
    dz = z .- zprev
  else
    mz = _reshape(mass_matrix * _vec(zprev), axes(zprev))
    dz = dt .* f(ztmp, p, tstep) .- mz
    z = zprev .+ dz
  end
  if has_destats(integrator)
    integrator.destats.nf += 1
  end

  # save residuals and proposed candidate 
  cache.dz = dz
  nlsolver.z = z

  nothing
end

@muladd function perform_step!(nlsolver::NLSolver{<:Union{NLFunctional,NLAnderson},true}, integrator, iter::Int)
  @unpack p,dt = integrator
  @unpack z,zprev,tmp,γ,cache = nlsolver
  @unpack dz,tstep,k = cache

  mass_matrix = integrator.f.mass_matrix
  f = nlsolve_f(integrator)

  # use z as temporary variable
  ztmp = z

  # evaluate function
  @.. ztmp = tmp + γ * zprev
  f(k, ztmp, p, tstep)
  if has_destats(integrator)
    integrator.destats.nf += 1
  end
  
  if mass_matrix == I
    @.. z = dt * k
    @.. dz = z - zprev
  else
    mul!(vec(ztmp), mass_matrix, vec(zprev))
    @.. dz = dt * k - ztmp
    @.. z = zprev + dz
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
  resize!(nlcache.gprev, i)
  resize!(nlcache.dzprev, i)

  # determine new maximum history
  max_history_old = length(Δgs)
  max_history = min(nlalg.max_history, nlalg.max_iter, i)

  resize!(nlcache.γs, max_history)
  resize!(nlcache.Δgs, max_history)
  if max_history > max_history_old
    for i in (max_history_old + 1):max_history
      nlcache.Δgs[i] = zero(nlcache.gprev)
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
  max_history_old = length(nlcache.Δgs)
  max_history = min(nlalg.max_history, nlalg.max_iter, i)

  resize!(nlcache.γs, max_history)
  resize!(nlcache.Δgs, max_history)

  if max_history != max_history_old
    nlcache.Q = typeof(nlcache.Q)(undef, i, max_history)
    nlcache.R = typeof(nlcache.R)(undef, max_history, max_history)
  end

  nothing
end