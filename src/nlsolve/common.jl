## build_nlsolver

function build_nlsolver(alg,nlalg::Union{NLFunctional,NLAnderson,NLNewton},u,rate_prototype,
                        uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,p,γ,c,
                        ::Val{true})
  # define unitless types
  uTolType = real(uBottomEltypeNoUnits)

  # define fields of non-linear solver
  z = similar(u); zprev = similar(u); tmp = similar(u)

  # build cache for non-linear solver
  dz = similar(u)
  tstep = zero(t)
  k = zero(rate_prototype)
  atmp = similar(u, uEltypeNoUnits)

  if nlalg isa NLNewton
    nf = nlsolve_f(f, alg)

    if islinear(f)
      du1 = rate_prototype
      uf = nothing
      jac_config = nothing
      linsolve = alg.linsolve(Val{:init},nf,u)
    else
      du1 = zero(rate_prototype)
      # if the algorithm specializes on split problems then use `nf`
      # we pass this `alg` here just for identification purpose, because get_uf would be overloaded in different repos
      uf = build_uf(alg,nf,t,p,Val(true))
      jac_config = build_jac_config(alg,nf,uf,du1,uprev,u,tmp,dz)
      linsolve = alg.linsolve(Val{:init},uf,u)
    end
    # TODO: check if the solver is iterative
    weight = similar(u)

    J, W = build_J_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits,Val(true))

    invγdt = inv(oneunit(dt) * one(uTolType))

    cache = NLNewtonCache(dz, tstep, k, atmp, J, W, true, dt, du1, uf, jac_config,
                          linsolve, weight, invγdt, (typeof(t))(nlalg.new_W_dt_cutoff))
  elseif nlalg isa NLFunctional
    cache = NLFunctionalCache(dz, tstep, k, atmp)
  elseif nlalg isa NLAnderson
    dzprev = similar(u)
    gprev = similar(u)

    max_history = min(nlalg.max_history, nlalg.maxiters, length(u))
    Δgs = [zero(u) for i in 1:max_history]
    Q = Matrix{uEltypeNoUnits}(undef, length(u), max_history)
    R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
    γs = Vector{uEltypeNoUnits}(undef, max_history)

    droptol = nlalg.droptol === nothing ? nothing : uEltypeNoUnits(nlalg.droptol)

    cache = NLAndersonCache(dz, tstep, k, atmp, gprev, dzprev, Δgs, Q, R, γs, 0, droptol)
  end

  # build non-linear solver
  η = one(uTolType)

  NLSolver{typeof(nlalg),true,typeof(u),uTolType,tTypeNoUnits,typeof(cache)}(
    z, zprev, tmp, uTolType(γ), tTypeNoUnits(c), nlalg, uTolType(nlalg.κ), η,
    uTolType(nlalg.fast_convergence_cutoff), nlalg.maxiters, 10_000, Convergence, cache)  
end

function build_nlsolver(alg,nlalg::Union{NLFunctional,NLAnderson,NLNewton},u,rate_prototype,
                        uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,p,γ,c,
                        ::Val{false})
  # define unitless types
  uTolType = real(uBottomEltypeNoUnits)

  # define fields of non-linear solver
  z = u; zprev = u; tmp = u

  # create cache of non-linear solver
  dz = u
  tstep = zero(t)

  if nlalg isa NLNewton
    nf = nlsolve_f(f, alg)
    # only use `nf` if the algorithm specializes on split eqs
    uf = build_uf(alg,nf,t,p,Val(false))

    J, W = build_J_W(alg,u,uprev,p,t,dt,f,uEltypeNoUnits,Val(false))

    invγdt = inv(oneunit(dt) * one(uTolType))

    cache = NLNewtonConstantCache(dz, tstep, J, W, true, dt, uf, invγdt,
                                  (typeof(t))(nlalg.new_W_dt_cutoff))
  elseif nlalg isa NLFunctional
    cache = NLFunctionalConstantCache(dz, tstep)
  elseif nlalg isa NLAnderson
    dzprev = u
    gprev = u
    
    max_history = min(nlalg.max_history, nlalg.maxiters, length(z))
    Δgs = Vector{typeof(u)}(undef, max_history)
    Q = Matrix{uEltypeNoUnits}(undef, length(u), max_history)
    R = Matrix{uEltypeNoUnits}(undef, max_history, max_history)
    γs = Vector{uEltypeNoUnits}(undef, max_history)

    droptol = nlalg.droptol === nothing ? nothing : uEltypeNoUnits(nlalg.droptol)

    cache = NLAndersonConstantCache(dz, tstep, gprev, dzprev, Δgs, Q, R, γs, 0, droptol)
  end

  # build non-linear solver
  η = one(uTolType)
  
  NLSolver{typeof(nlalg),false,typeof(u),uTolType,tTypeNoUnits,typeof(cache)}(
    z, zprev, tmp, uTolType(γ), tTypeNoUnits(c), nlalg, uTolType(nlalg.κ), η,
    uTolType(nlalg.fast_convergence_cutoff), nlalg.maxiters, 10_000, Convergence, cache)
end

## norm_of_residuals

function norm_of_residuals(nlsolver::NLSolver{<:Union{NLFunctional,NLAnderson,NLNewton},true},
                           integrator)
  @unpack t,opts = integrator
  @unpack z,zprev,cache = nlsolver
  @unpack dz,atmp = cache

  calculate_residuals!(atmp, dz, zprev, z, opts.abstol, opts.reltol, opts.internalnorm, t)
  opts.internalnorm(atmp, t)
end

function norm_of_residuals(nlsolver::NLSolver{<:Union{NLFunctional,NLAnderson,NLNewton},false},
                           integrator)
  @unpack t,opts = integrator
  @unpack z,zprev,cache = nlsolver
  @unpack dz = cache

  atmp = calculate_residuals(dz, zprev, z, opts.abstol, opts.reltol, opts.internalnorm, t)
  opts.internalnorm(atmp, t)
end