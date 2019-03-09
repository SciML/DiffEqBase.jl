mutable struct DEStats
  nf::Int
  nf2::Int
  nw::Int
  nsolve::Int
  njacs::Int
  nnonliniter::Int
  nnonlinconvfail::Int
  ncondition::Int
  naccept::Int
  nreject::Int
  maxeig::Float64
end

DEStats(x=-1) = DEStats(ntuple(_->x, 10)..., 0.0)

function Base.show(io::IO, s::DEStats)
  println(io, summary(s))
  println(io, "Number of function 1 evaluations: $(s.nf)")
  println(io, "Number of function 2 evaluations: $(s.nf2)")
  println(io, "Number of W matrix evaluations: $(s.nw)")
  println(io, "Number of linear solves:        $(s.nsolve)")
  println(io, "Number of Jacobians created:        $(s.njacs)")
  println(io, "Number of nonlinear solver iterations:        $(s.nnonliniter)")
  println(io, "Number of nonlinear solver convergence failures:        $(s.nnonlinconvfail)")
  println(io, "Number of rootfind condition calls:        $(s.ncondition)")
  println(io, "Number of accepted steps:       $(s.naccept)")
  print(io,   "Number of rejected steps:       $(s.nreject)")
  iszero(s.maxeig) || print(io, "\nMaximum eigenvalue recorded:    $(s.maxeig)")
end
