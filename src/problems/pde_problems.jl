"""
$(TYPEDEF)

TODO
"""
struct PDEProblem{P,E,S} <: AbstractPDEProblem
  prob::P
  extrapolation::E
  space::S
end
