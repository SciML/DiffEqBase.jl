using DiffEqBase
mutable struct DummySolution
  retcode
end

mutable struct DummyIntegrator{Alg, IIP, U, T} <: DiffEqBase.DEIntegrator{Alg, IIP, U, T}
  t
  dt
  tdir
  tstops
  sol

  DummyIntegrator() = new{Bool,Bool,Bool,Bool}(0,1,1,[],DummySolution(:Default))
end

function DiffEqBase.add_tstop!(integrator::DummyIntegrator,t)
  integrator.tdir * (t - integrator.t) < 0 && error("Tried to add a tstop that is behind the current time. This is strictly forbidden")
  push!(integrator.tstops,t)
end

function DiffEqBase.step!(integrator::DummyIntegrator)
  t_next = integrator.t + integrator.dt
  if ! isempty(integrator.tstops) && integrator.tstops[1] < t_next
    t_next = integrator.tstops[1]
  end
  integrator.t = t_next
end

function step_dt!(integrator, args...)
  t = integrator.t
  step!(integrator, args...)
  integrator.t - t
end

integrator = DummyIntegrator()
@test step_dt!(integrator, 1.5) == 2
@test step_dt!(integrator, 1.5, true) == 1.5
@test_throws ErrorException step!(integrator, -1)
