struct DEStat{E}
  nf::UInt
  njac::UInt
  nsolve::UInt
  naccept::UInt
  nreject::UInt
  maxeig::E
end

DEStat(eig) = DEStat(ntuple(_->UInt(0), 5)..., eig)

function Base.show(io::IO, s::DEStat)
  println(io, summary(s))
  printstyled(io, "Number of function evaluations:  $(s.nf)")
  printstyled(io, "Number of Jacobian evaluations:  $(s.njac)")
  printstyled(io, "Number of linear solves:         $(s.nsolve)")
  printstyled(io, "Number of accepted steps:        $(s.accepted)")
  printstyled(io, "Number of rejected steps:        $(s.rejected)")
  printstyled(io, "Maximum eigenvalue recorded:     $(s.maxeig)")
end
