type Callback{F1,F2} <: DECallback
  condition::F1
  affect!::F2
  rootfind::Bool
  interp_points::Int
  save_positions::Tuple{Bool,Bool}
end
