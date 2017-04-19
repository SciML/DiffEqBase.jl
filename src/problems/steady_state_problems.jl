# Mu' = f
type SteadyStateProblem{uType,isinplace,F,MM} <: AbstractSteadyStateProblem{uType,isinplace}
  f::F
  u0::uType
  mass_matrix::MM
end

function SteadyStateProblem(f,u0; iip = isinplace(f,3),mass_matrix=I)
  SteadyStateProblem{typeof(u0),iip,typeof(f),typeof(mass_matrix)}(f,u0,mass_matrix)
end
