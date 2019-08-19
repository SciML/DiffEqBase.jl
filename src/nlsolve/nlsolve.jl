"""
    nlsolve!(nlsolver::AbstractNLSolver, integrator)

Solve
```math
G(z) = dt⋅f(tmp + γ⋅z, p, t + c⋅h) - z = 0
```
where `dt` is the step size and `γ` and `c` are constants, and return the solution `z`.
"""
function nlsolve!(nlsolver::AbstractNLSolver, integrator)
  preamble!(nlsolver, integrator)

  while get_status(nlsolver) === SlowConvergence
    # (possibly modify and) accept step
    loopheader!(nlsolver, integrator)

    # compute next iterate
    perform_step!(nlsolver, integrator)
  
    # check convergence and divergence criteria
    loopfooter!(nlsolver, integrator)
  end

  postamble!(nlsolver, integrator)
end

## default implementations for NLSolver

function preamble!(nlsolver::NLSolver, integrator)
  nlsolver.iter = 0
  if nlsolver.maxiters == 0
    nlsolver.status = MaxIterReached
    return
  end
  
  nlsolver.status = SlowConvergence
  nlsolver.η = initial_η(nlsolver, integrator)

  initialize_chache!(nlsolver.cache, nlsolver, integrator)

  nothing
end

initial_η(nlsolver::NLSolver, integrator) = nlsolver.ηold

function loopheader!(nlsolver::NLSolver, integrator)
  # apply step
  # adjust next iterate
  # can be used, e.g., for Anderson acceleration
  adjust_step!(nlsolver, integrator)
  
  # apply step
  apply_step!(nlsolver, integrator)

  # update statistics
  nlsolver.iter += 1
  if has_destats(integrator)
    integrator.destats.nnonliniter += 1
  end

  nothing
end

function apply_step!(nlsolver::NLSolver{algType,iip}, integrator) where {algType,iip}
  if iip
    recursivecopy!(nlsolver.zprev, nlsolver.z)
  else
    nlsolver.zprev = nlsolver.z
  end

  nothing
end

loopfooter!(nlsolver::NLSolver, integrator) = check_status!(nlsolver, integrator)

function check_status!(nlsolver::NLSolver, integrator)
  nlsolver.status = check_status(nlsolver, integrator)
  nothing
end

function check_status(nlsolver::NLSolver, integrator)
  @unpack iter,maxiters,κ,fast_convergence_cutoff = nlsolver

  # compute norm of residuals and cache previous value
  iter > 1 && (ndzprev = ndz)
  ndz = norm_of_residuals(nlsolver, integrator)

  # check for convergence
  if iter > 1
    Θ = ndz / ndzprev
    η = Θ / (1 - Θ)
    nlsolver.η = η
  else
    η = nlsolver.η
  end
  if iszero(ndz) || (η * ndz < κ && (iter > 1 || !iszero(integrator.success_iter)))
    if η < nlsolver.fast_convergence_cutoff
      return FastConvergence
    else
      return Convergence
    end
  end

  # check for divergence (not in initial step)
  if iter > 1
    # divergence
    if Θ > 1
      return Divergence
    end

    # very slow convergence
    if ndz * Θ^(maxiters - iter) > κ * (1 - Θ)
      return VerySlowConvergence
    end
  end

  # check number of iterations
  if iter >= maxiters
    return MaxIterReached
  end

  SlowConvergence
end  

function norm_of_residuals(nlsolver::NLSolver, integrator)
  @unpack t,opts = integrator
  @unpack z,zprev = nlsolver

  atmp = calculate_residuals(zprev, z, opts.abstol, opts.reltol, opts.internalnorm, t)
  opts.internalnorm(atmp, t)
end

function postamble!(nlsolver::NLSolver, integrator)
  fail_convergence = nlsolvefail(nlsolver)
  if fail_convergence && has_destats(integrator)
      integrator.destats.nnonlinconvfail += 1
  end
  integrator.force_stepfail = fail_convergence

  # return previous iterate
  # (we do not have any convergence guarantees on the subsequent one)
  nlsolver.zprev
end
