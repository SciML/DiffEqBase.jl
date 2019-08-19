"""
    nlsolve!(nlsolver::NLSolver, integrator)

Solve
```math
G(z) = dt⋅f(tmp + γ⋅z, p, t + c⋅h) - z = 0
```
where `dt` is the step size and `γ` and `c` are constants, and return the solution `z`.
"""
function nlsolve!(nlsolver::NLSolver, integrator)
  @unpack max_iter,κ,fast_convergence_cutoff = nlsolver

  preamble!(nlsolver, integrator)

  local ndz
  η = initial_η(nlsolver, integrator)
  iter = 0
  while true
    # check number of iterations
    if iter >= max_iter
      nlsolver.status = MaxIterReached
      break
    end

    # apply step
    loopheader!(nlsolver, integrator, iter)

    # compute next iterate
    iter += 1
    if has_destats(integrator)
      integrator.destats.nnonliniter += 1
    end
    perform_step!(nlsolver, integrator, iter)
  
    # compute norm of residuals and cache previous value
    iter > 1 && (ndzprev = ndz)
    ndz = norm_of_residuals(nlsolver, integrator)

    # check for convergence
    if iter > 1
      Θ = ndz / ndzprev
      η = Θ / (1 - Θ)
    end
    if iszero(ndz) || (η * ndz < κ && (iter > 1 || !iszero(integrator.success_iter)))
      nlsolver.status = η < fast_convergence_cutoff ? FastConvergence : Convergence
      break
    end

    # check for divergence (not in initial step)
    if iter > 1
      # divergence
      if Θ > 1
        nlsolver.status = Divergence
        break
      end

      # very slow convergence
      if ndz * Θ^(max_iter - iter) > κ * (1 - Θ)
        nlsolver.status = VerySlowConvergence
        break
      end
    end
  end

  fail_convergence = nlsolvefail(nlsolver)
  if fail_convergence && has_destats(integrator)
      integrator.destats.nnonlinconvfail += 1
  end
  integrator.force_stepfail = fail_convergence

  nlsolver.ηold = η
  nlsolver.nl_iters = iter
    
  # return previous iterate
  # (we do not have any convergence guarantees on the subsequent one)
  nlsolver.zprev
end

## default implementations

initial_η(nlsolver::NLSolver, integrator) = nlsolver.ηold

preamble!(nlsolver::NLSolver, integrator) = nothing

loopheader!(nlsolver::NLSolver, integrator, iter::Int) = apply_step!(nlsolver, integrator)

function apply_step!(nlsolver::NLSolver{algType,iip}, integrator) where {algType,iip}
  if iip
    recursivecopy!(nlsolver.zprev, nlsolver.z)
  else
    nlsolver.zprev = nlsolver.z
  end

  nothing
end

function norm_of_residuals(nlsolver::NLSolver, integrator)
  @unpack t,opts = integrator
  @unpack z,zprev = nlsolver

  atmp = calculate_residuals(zprev, z, opts.abstol, opts.reltol, opts.internalnorm, t)
  opts.internalnorm(atmp, t)
end