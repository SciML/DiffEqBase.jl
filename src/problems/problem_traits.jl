is_diagonal_noise{uType,tType,isinplace,ND}(prob::AbstractRODEProblem{uType,tType,isinplace,ND}) = ND <: Void

isinplace{uType,tType,iip}(prob::AbstractODEProblem{uType,tType,iip}) = iip
isinplace{uType,iip}(prob::AbstractSteadyStateProblem{uType,iip}) = iip
