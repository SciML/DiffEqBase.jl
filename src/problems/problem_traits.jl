"""
    is_diagonal_noise(prob::DEProblem)

TODO
"""
is_diagonal_noise(prob::DEProblem) = false
is_diagonal_noise(prob::AbstractRODEProblem{uType,tType,iip,Nothing}) where {uType,tType,iip} = true
is_diagonal_noise(prob::AbstractSDDEProblem{uType,tType,lType,iip,Nothing}) where {uType,tType,lType,iip,ND} = true

"""
    isinplace(prob::DEProblem)

Determine whether the function of the given problem operates in place or not.
"""
function isinplace(prob::DEProblem) end
isinplace(prob::AbstractLinearProblem{bType,iip}) where {bType,iip} = iip
isinplace(prob::AbstractNonlinearProblem{uType,iip}) where {uType,iip} = iip
isinplace(prob::AbstractQuadratureProblem{iip}) where {iip} = iip
isinplace(prob::AbstractODEProblem{uType,tType,iip}) where {uType,tType,iip} = iip
isinplace(prob::AbstractSteadyStateProblem{uType,iip}) where {uType,iip} = iip
isinplace(prob::AbstractRODEProblem{uType,tType,iip,ND}) where {uType,tType,iip,ND} = iip
isinplace(prob::AbstractDDEProblem{uType,tType,lType,iip}) where {uType,tType,lType,iip} = iip
isinplace(prob::AbstractDAEProblem{uType,duType,tType,iip}) where {uType,duType,tType,iip} = iip
isinplace(prob::AbstractNoiseProblem) = isinplace(prob.noise)
isinplace(::SplitFunction{iip}) where iip = iip
isinplace(prob::AbstractSDDEProblem{uType,tType,lType,iip,ND}) where {uType,tType,lType,iip,ND} = iip
