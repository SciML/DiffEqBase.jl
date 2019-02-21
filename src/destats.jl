mutable struct DEStat{E}
  nf::Int
  nw::Int
  nsolve::Int
  naccept::Int
  nreject::Int
  maxeig::E
end

DEStat(eig) = DEStat(ntuple(_->0, 5)..., eig)

function Base.show(io::IO, s::DEStat)
  println(io, summary(s))
  println(io, "Number of function evaluations: $(s.nf)")
  println(io, "Number of W matrix evaluations: $(s.nw)")
  println(io, "Number of linear solves:        $(s.nsolve)")
  println(io, "Number of accepted steps:       $(s.naccept)")
  print(io,   "Number of rejected steps:       $(s.nreject)")
  s.maxeig !== nothing && print(io, "\nMaximum eigenvalue recorded:    $(s.maxeig)")
end
