struct PDEProblem{P,E,S} <: DiffEqBase.DEProblem
  prob::P
  extrapolation::E
  space::S
end
