function problem_new_parameters(prob::ODEProblem,p;kwargs...)
  uEltype = eltype(p)
  u0 = [uEltype(prob.u0[i]) for i in 1:length(prob.u0)]
  tspan = (uEltype(prob.tspan[1]),uEltype(prob.tspan[2]))
  ODEProblem{isinplace(prob)}(prob.f,u0,tspan,p,prob.problem_type;
  callback = prob.callback,
  kwargs...)
end

function problem_new_parameters(prob::BVProblem,p;kwargs...)
  uEltype = eltype(p)
  u0 = [uEltype(prob.u0[i]) for i in 1:length(prob.u0)]
  tspan = (uEltype(prob.tspan[1]),uEltype(prob.tspan[2]))
  BVProblem{isinplace(prob)}(prob.f,prob.bc,u0,tspan,p,prob.problem_type;
  callback = prob.callback,
  kwargs...)
end

function problem_new_parameters(prob::DAEProblem,p;kwargs...)
  uEltype = eltype(p)
  u0 = [uEltype(prob.u0[i]) for i in 1:length(prob.u0)]
  du0 = [uEltype(prob.du0[i]) for i in 1:length(prob.du0)]
  tspan = (uEltype(prob.tspan[1]),uEltype(prob.tspan[2]))
  DAEProblem{isinplace(prob)}(prob.f,du0,u0,tspan,p;
  differential_vars=prob.differential_vars,
  callback = prob.callback,
  kwargs...)
end

function problem_new_parameters(prob::DDEProblem,p;kwargs...)
  uEltype = eltype(p)
  u0 = [uEltype(prob.u0[i]) for i in 1:length(prob.u0)]
  tspan = (uEltype(prob.tspan[1]),uEltype(prob.tspan[2]))
  DDEProblem{isinplace(prob)}(prob.f,u0,prob.h,tspan,p;
                              constant_lags = prob.constant_lags,
                              dependent_lags = prob.dependent_lags,
  callback = prob.callback,
  neutral = prob.neutral,
  kwargs...)
end

function problem_new_parameters(prob::SDEProblem,p;kwargs...)
  uEltype = eltype(p)
  u0 = [uEltype(prob.u0[i]) for i in 1:length(prob.u0)]
  tspan = (uEltype(prob.tspan[1]),uEltype(prob.tspan[2]))
  SDEProblem{isinplace(prob)}(prob.f,prob.g,u0,tspan,p;
  noise_rate_prototype = prob.noise_rate_prototype,
  noise= prob.noise, seed = prob.seed,
  callback = prob.callback,
  kwargs...)
end

function problem_new_parameters(prob::MonteCarloProblem,p)
  MonteCarloProblem(problem_new_parameters(prob.prob,p),
                    prob.prob_func,prob.output_func,prob.reduction,prob.u_init)
end
