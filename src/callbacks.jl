immutable Callback{F1,F2,T,T2} <: DECallback
  condition::F1
  affect!::F2
  rootfind::Bool
  interp_points::Int
  save_positions::Tuple{Bool,Bool}
  abstol::T
  reltol::T2
end

Callback(condition,affect!,
         rootfind,interp_points,
         save_positions;
         abstol=0,reltol=1e-10) = Callback(condition,affect!,
                    rootfind,interp_points,
                    save_positions,abstol,reltol)
