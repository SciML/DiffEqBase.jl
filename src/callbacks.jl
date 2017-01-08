immutable Callback{F1,F2,F3,T,T2} <: DECallback
  condition::F1
  affect!::F2
  affect_neg!::F3
  rootfind::Bool
  interp_points::Int
  save_positions::Tuple{Bool,Bool}
  abstol::T
  reltol::T2
end

Callback(condition,affect!,
         rootfind,
         save_positions;
         affect_neg! = affect!,
         interp_points=10,
         abstol=1e-14,reltol=0) = Callback(condition,affect!,affect_neg!,
                    rootfind,interp_points,
                    save_positions,abstol,reltol)
