is_diagonal_noise{uType,tType,isinplace,ND}(prob::AbstractRODEProblem{uType,tType,isinplace,ND}) = ND <: Void

isinplace{uType,tType,iip}(prob::AbstractODEProblem{uType,tType,iip}) = iip
isinplace{uType,iip}(prob::AbstractSteadyStateProblem{uType,iip}) = iip
isinplace{uType,tType,iip,ND}(prob::AbstractRODEProblem{uType,tType,iip,ND}) = iip
isinplace{uType,tType,lType,iip}(prob::AbstractDDEProblem{uType,tType,lType,iip}) = iip
isinplace{uType,duType,tType,iip}(prob::AbstractDAEProblem{uType,duType,tType,iip}) = iip
