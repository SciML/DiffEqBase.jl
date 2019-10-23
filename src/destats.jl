"""
$(TYPEDEF)

TODO
"""
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

DEStats(x::Int = -1) = DEStats(x, x, x, x, x, x, x, x, x, x, 0.0)

function Base.show(io::IO, s::DEStats)
  println(io, summary(s))
  @printf io "%-50s %-d\n" "Number of function 1 evaluations:" s.nf
  @printf io "%-50s %-d\n" "Number of function 2 evaluations:" s.nf2
  @printf io "%-50s %-d\n" "Number of W matrix evaluations:" s.nw
  @printf io "%-50s %-d\n" "Number of linear solves:" s.nsolve
  @printf io "%-50s %-d\n" "Number of Jacobians created:" s.njacs
  @printf io "%-50s %-d\n" "Number of nonlinear solver iterations:" s.nnonliniter
  @printf io "%-50s %-d\n" "Number of nonlinear solver convergence failures:" s.nnonlinconvfail
  @printf io "%-50s %-d\n" "Number of rootfind condition calls:" s.ncondition
  @printf io "%-50s %-d\n" "Number of accepted steps:"  s.naccept
  @printf io "%-50s %-d" "Number of rejected steps:" s.nreject
  iszero(s.maxeig) || @printf io "\n%-50s %-d" "Maximum eigenvalue recorded:" s.maxeig
end
