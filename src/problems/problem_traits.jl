is_diagonal_noise{uType,tType,isinplace,NoiseClass,ND}(prob::AbstractRODEProblem{uType,tType,isinplace,NoiseClass,ND}) = ND <: Void

isinplace{uType,tType,iip}(prob::AbstractODEProblem{uType,tType,iip}) = iip
