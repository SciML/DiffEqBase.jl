is_diagonal_noise{uType,tType,isinplace,ND}(prob::AbstractRODEProblem{uType,tType,isinplace,ND}) = ND <: Void

isinplace{uType,tType,iip}(prob::AbstractODEProblem{uType,tType,iip}) = iip
isinplace{uType,iip}(prob::AbstractSteadyStateProblem{uType,iip}) = iip
isinplace{uType,tType,iip,ND}(prob::AbstractRODEProblem{uType,tType,iip,ND}) = iip
isinplace{uType,tType,lType,iip}(prob::AbstractDDEProblem{uType,tType,lType,iip}) = iip
isinplace{uType,duType,tType,iip}(prob::AbstractDAEProblem{uType,duType,tType,iip}) = iip
isinplace(prob::AbstractNoiseProblem) = isinplace(prob.noise)

### Displays

Juno.@render Juno.Inline x::DEProblem begin
  fields = fieldnames(typeof(x))
  Juno.LazyTree(typeof(x), () -> [Juno.SubTree(Juno.Text("$f → "), Juno.getfield′(x, f)) for f in fields])
end

Base.summary(prob::DEProblem) = string(parameterless_type(prob)," with uType ",typeof(prob.u0)," and tType ",typeof(prob.tspan[1]),". In-place: ",isinplace(prob))
Base.summary{uType,iip}(prob::AbstractSteadyStateProblem{uType,iip}) = string(parameterless_type(prob)," with uType ",uType)
Base.summary(prob::AbstractNoiseProblem) = string(parameterless_type(prob)," with WType ",typeof(prob.noise.W[1])," and tType ",typeof(prob.tspan[1]),". In-place: ",isinplace(prob))
function Base.show(io::IO, A::DEProblem)
  println(io,summary(A))
  print(io,"timespan: ")
  show(io,A.tspan)
  println(io)
  print(io,"u0: ")
  show(io, A.u0)
end
function Base.show(io::IO, A::AbstractNoiseProblem)
  println(io,summary(A))
  print(io,"timespan: ")
  show(io,A.tspan)
  println(io)
end
function Base.show(io::IO, A::AbstractDAEProblem)
  println(io,summary(A))
  print(io,"timespan: ")
  show(io,A.tspan)
  println(io)
  print(io,"u0: ")
  show(io, A.u0)
  println(io)
  print(io,"du0: ")
  show(io, A.du0)
end
function Base.show(io::IO, A::AbstractSteadyStateProblem)
  println(io,summary(A))
  print(io,"u0: ")
  show(io, A.u0)
end
function Base.display(io::IO, A::DEProblem)
  println(io,summary(A))
  print(io,"timespan: ")
  display(io,A.tspan)
  println(io)
  print(io,"u0: ")
  display(io, A.u0)
end
function Base.display(io::IO, A::AbstractNoiseProblem)
  println(io,summary(A))
  print(io,"timespan: ")
  display(io,A.tspan)
  println(io)
end
function Base.display(io::IO, A::AbstractDAEProblem)
  println(io,summary(A))
  print(io,"timespan: ")
  display(io,A.tspan)
  println(io)
  print(io,"u0: ")
  display(io, A.u0)
  println(io)
  print(io,"du0: ")
  display(io, A.du0)
end
function Base.print(io::IO,A::DEProblem)
  show(io,A)
end
function Base.println(io::IO,A::DEProblem)
  show(io,A)
end
Base.print(A::DEProblem) = print(STDOUT,A)
Base.println(A::DEProblem) = println(STDOUT,A)

Base.summary(prob::AbstractMonteCarloProblem) = string(DiffEqBase.parameterless_type(prob)," with problem ",DiffEqBase.parameterless_type(prob.prob))
function Base.show(io::IO, A::AbstractMonteCarloProblem)
  println(io,summary(A))
end
function Base.display(io::IO, A::AbstractMonteCarloProblem)
  println(io,summary(A))
end
