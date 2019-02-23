mutable struct DEStats
  nf::Int
  nw::Int
  nsolve::Int
  naccept::Int
  nreject::Int
  maxeig::Float64
end

DEStats() = DEStats(ntuple(_->0, 5)..., 0.0)

function Base.show(io::IO, s::DEStats)
  println(io, summary(s))
  println(io, "Number of function evaluations: $(s.nf)")
  println(io, "Number of W matrix evaluations: $(s.nw)")
  println(io, "Number of linear solves:        $(s.nsolve)")
  println(io, "Number of accepted steps:       $(s.naccept)")
  print(io,   "Number of rejected steps:       $(s.nreject)")
  iszero(s.maxeig) || print(io, "\nMaximum eigenvalue recorded:    $(s.maxeig)")
end
