using DiffEqBase
mutable struct DummySolution
  retcode
end

mutable struct DummyIntegrator{Alg, IIP, U, T} <: DiffEqBase.DEIntegrator{Alg, IIP, U, T}
  uprev::U
  tprev::T
  u::U
  t::T
  dt::T
  tdir
  tstops
  sol::DummySolution

  DummyIntegrator() = new{Bool,Bool,Vector{Float64},Float64}([0.0],0,[0.0],0,1,1,[],DummySolution(:Default))
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
  integrator.uprev .= integrator.u
  integrator.u[1] += 2 * (t_next - integrator.t)
  integrator.tprev = integrator.t
  integrator.t = t_next
end

function step_dt!(integrator, args...)
  t = integrator.t
  step!(integrator, args...)
  integrator.t - t
end

function DiffEqBase.done(integrator::DummyIntegrator)
  integrator.t > 10
end

integrator = DummyIntegrator()
@test step_dt!(integrator, 1.5) == 2
@test step_dt!(integrator, 1.5, true) == 1.5
@test_throws ErrorException step!(integrator, -1)

for (u, t) in tuples(DummyIntegrator())
  @test u[1] == 2*t
end
@test eltype(collect(tuples(DummyIntegrator()))) == Tuple{Vector{Float64},Float64}

for (uprev, tprev, u, t) in intervals(DummyIntegrator())
  @test u[1] - uprev[1] == 2
  @test t - tprev == 1
end
@test eltype(collect(intervals(DummyIntegrator()))) == Tuple{Vector{Float64},Float64,Vector{Float64},Float64}