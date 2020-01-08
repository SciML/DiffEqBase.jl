"""
    promote_tspan(tspan)

Convert the `tspan` field of a `DEProblem` to a `(tmin, tmax)` tuple, where both
elements are of the same type. If `tspan` is a function, returns it as-is.
"""
promote_tspan((t1,t2)::Tuple{T,S}) where {T,S} = promote(t1, t2)
promote_tspan(tspan::Number) = (zero(tspan),tspan)
promote_tspan(tspan::Nothing) = (nothing,nothing)
promote_tspan(tspan::Function) = tspan

### Displays

Base.summary(prob::DEProblem) = string(TYPE_COLOR, nameof(typeof(prob)),
                                       NO_COLOR, " with uType ",
                                       TYPE_COLOR, typeof(prob.u0),
                                       NO_COLOR, " and tType ",
                                       TYPE_COLOR,
                                       typeof(prob.tspan) <: Function ?
                                       "Unknown" : (prob.tspan === nothing ?
                                       "Nothing" : typeof(prob.tspan[1])),
                                       NO_COLOR, ". In-place: ",
                                       TYPE_COLOR, isinplace(prob),
                                       NO_COLOR)


Base.summary(prob::AbstractLinearProblem) = string(TYPE_COLOR, nameof(typeof(prob)),
                                                   NO_COLOR, ". In-place: ",
                                                   TYPE_COLOR, isinplace(prob),
                                                   NO_COLOR)
function Base.show(io::IO, A::AbstractLinearProblem)
  show(io,summary(A))
  print(io,"b: ")
  show(io, A.b)
end

Base.summary(prob::AbstractNonlinearProblem) = string(
                                       TYPE_COLOR, nameof(typeof(prob)),
                                       NO_COLOR, ". In-place: ",
                                       TYPE_COLOR, isinplace(prob),
                                       NO_COLOR)
function Base.show(io::IO, A::AbstractNonlinearProblem)
  println(io,summary(A))
  print(io,"u0: ")
  show(io, A.u0)
end

Base.summary(prob::AbstractQuadratureProblem) = string(
                                                       TYPE_COLOR, nameof(typeof(prob)),
                                                       NO_COLOR, ". In-place: ",
                                                       TYPE_COLOR, isinplace(prob),
                                                       NO_COLOR)
function Base.show(io::IO, A::AbstractQuadratureProblem)
  println(io,summary(A))
end

Base.summary(prob::AbstractSteadyStateProblem{uType,iip}) where {uType,iip} = string(nameof(typeof(prob))," with uType ",uType)
Base.summary(prob::AbstractNoiseProblem) = string(nameof(typeof(prob))," with WType ",typeof(prob.noise.W[1])," and tType ",typeof(prob.tspan[1]),". In-place: ",isinplace(prob))
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

Base.summary(prob::AbstractEnsembleProblem) = string(
nameof(typeof(prob))," with problem ",nameof(typeof(prob.prob)))
Base.show(io::IO, A::AbstractEnsembleProblem) = print(io,summary(A))
TreeViews.hastreeview(x::DiffEqBase.DEProblem) = true
function TreeViews.treelabel(io::IO,x::DiffEqBase.DEProblem,
                             mime::MIME"text/plain" = MIME"text/plain"())
  show(io,mime,Base.Text(Base.summary(x)))
end

struct NullParameters end
Base.getindex(::NullParameters,i...) = error("Parameters were indexed but the parameters are `nothing`. You likely forgot to pass in parameters to the DEProblem!")

function Base.show(io::IO, A::AbstractPDEProblem)
  println(io,summary(A.prob))
  println(io)
end
Base.summary(prob::AbstractPDEProblem) = string(DiffEqBase.TYPE_COLOR,
                                                nameof(typeof(prob)),
                                                DiffEqBase.NO_COLOR)
