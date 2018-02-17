mutable struct DummyIntegrator <: DEIntegrator
  t
  dt
  tstops

  DummyIntegrator() = new(0,1,[])
end

function DiffEqBase.add_tstop!(integrator::DummyIntegrator,t)
  integrator.dt * (t - integrator.t) < 0 && error("Tried to add a tstop that is behind the current time. This is strictly forbidden")
  push!(integrator.tstops,t)
end

function DiffEqBase.step!(integrator::DummyIntegrator)
  t_next = integrator.t + integrator.dt
  if ! isempty(integrator.tstops) && integrator.tstops[1] < t_next
    t_next = integrator.tstops[1]
  end
  integrator.t = t_next
end

integrator = DummyIntegrator()
@test step!(integrator, 1.5) == 2
@test step!(integrator, 1.5, true) == 1.5
