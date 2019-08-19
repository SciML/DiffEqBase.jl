## initial_η

initial_η(nlsolver::NLSolver{NLNewton}, integrator) =
  max(nlsolver.ηold, eps(eltype(integrator.opts.reltol)))^(0.8)

## preamble!

@muladd function initialize_cache!(nlcache::NLNewtonConstantCache,
                                   nlsolver::NLSolver{<:NLNewton,false}, integrator)
  @unpack dt = integrator

  nlcache.invγdt = inv(dt * nlsolver.γ)
  nlcache.tstep = integrator.t + nlsolver.c * dt 

  nothing
end

@muladd function initialize_cache!(nlcache::NLNewtonCache,
                                   nlsolver::NLSolver{<:NLNewton,true}, integrator)
  @unpack u,uprev,t,dt,opts = integrator
  @unpack weight = nlcache

  nlcache.invγdt = inv(dt * nlsolver.γ)
  nlcache.tstep = integrator.t + nlsolver.c * dt 
  calculate_residuals!(weight, fill!(weight, one(eltype(u))), uprev, u,
                       opts.abstol, opts.reltol, opts.internalnorm, t)
  
  nothing
end

## perform_step!

"""
    perform_step!(nlsolver::NLSolver{<:NLNewton}, integrator)

Compute next iterate of numerically stable modified Newton iteration
that is specialized for implicit methods (see [^HS96] and [^HW96]).

It solves
```math
G(z) = dt⋅f(tmp + γ⋅z, p, t + c⋅dt) - z = 0
```
by iterating
```math
(I + (dt⋅γ)J) Δᵏ = dt*f(tmp + γ⋅zᵏ, p, t + c⋅dt) - zᵏ
zᵏ⁺¹ = zᵏ + Δᵏ
```
or, by utilizing a transformation,
```math
W Δᵏ = f(tmp + γ⋅zᵏ, p, t + c⋅dt)/γ - zᵏ/(dt⋅γ)
zᵏ⁺¹ = zᵏ + Δᵏ/(dt⋅γ)
```
where `W = M/(dt⋅γ) - J`, `M` is the mass matrix, `dt` is the step size, `γ` and
`c` are constants, `J` is the Jacobian matrix. This transformation occurs since `c*J` is
O(n^2), while `c*M` is usually much sparser. In the most common case, `M=I`, we
have that `c*M` is O(1) for `I isa UniformScaling`.

[^HS96]: M.E.Hoseaa and L.F.Shampine, "Analysis and implementation of TR-BDF2",
Applied Numerical Mathematics, Volume 20, Issues 1–2, February 1996, Pages
21-37.
[doi:10.1016/0168-9274(95)00115-8](https://doi.org/10.1016/0168-9274(95)00115-8)

[^HW96]: Ernst Hairer and Gerhard Wanner, "Solving Ordinary Differential
Equations II, Springer Series in Computational Mathematics. ISBN
978-3-642-05221-7. Section IV.8.
[doi:10.1007/978-3-642-05221-7](https://doi.org/10.1007/978-3-642-05221-7)
"""
@muladd function perform_step!(nlsolver::NLSolver{<:NLNewton,false}, integrator)
  @unpack p,dt = integrator
  @unpack zprev,tmp,γ,cache = nlsolver
  @unpack tstep,W,invγdt = cache

  # precalculations
  mass_matrix = integrator.f.mass_matrix
  f = nlsolve_f(integrator)

  # evaluate function
  ztmp = @.. tmp + γ * zprev
  if mass_matrix === I
    ztmp = (dt .* f(ztmp, p, tstep) .- zprev) .* invγdt
  else
    ztmp = (dt .* f(ztmp, p, tstep) .- mass_matrix * zprev) .* invγdt
  end
  if has_destats(integrator)
    integrator.destats.nf += 1
  end

  dz = _reshape(W \ _vec(ztmp), axes(ztmp))
  if has_destats(integrator)
    integrator.destats.nsolve += 1
  end

  cache.dz = dz
  nlsolver.z = @.. zprev - dz

  nothing
end

@muladd function perform_step!(nlsolver::NLSolver{<:NLNewton,true}, integrator)
  @unpack p,dt = integrator
  @unpack z,zprev,tmp,γ,iter,cache = nlsolver
  @unpack dz,tstep,k,W,new_W,linsolve,weight,invγdt = cache
  
  mass_matrix = integrator.f.mass_matrix
  f = nlsolve_f(integrator)

  # use z as temporary variable
  ztmp = z

  # evaluate function
  @.. ztmp = tmp + γ * zprev
  if W isa AbstractDiffEqLinearOperator
    update_coefficients!(W, ztmp, p, tstep)
  end

  f(k, ztmp, p, tstep)
  if has_destats(integrator)
    integrator.destats.nf += 1
  end

  if mass_matrix === I
    @.. ztmp = (dt * k - zprev) * invγdt
  else
    mul!(vec(z), mass_matrix, vec(zprev))
    @.. ztmp = (dt * k - z) * invγdt
  end

  linsolve(vec(dz), W, vec(ztmp), iter == 1 && new_W; Pl=ScaleVector(weight, true), Pr=ScaleVector(weight, false), tol=integrator.opts.reltol)
  if has_destats(integrator)
    integrator.destats.nsolve += 1
  end

  # compute next iterate
  @.. z = zprev - dz

  nothing
end

## traits

isnewton(::NLSolver{<:NLNewton}) = true

## resize!

function Base.resize!(nlcache::NLNewtonCache, nlsolver, integrator, i::Int)
  resize!(nlcache.dz, i)
  resize!(nlcache.k, i)
  resize!(nlcache.atmp, i)
  resize!(nlcache.du1, i)
  if jac_config !== nothing
    nlcache.jac_config = resize_jac_config!(nlcache.jac_config, i)
  end
  resize!(nlcache.weight, i)

  # have to be implemented for integrator
  resize_J!(nlcache, integrator, i)
  resize_W!(nlcache, integrator, i)

  nothing
end