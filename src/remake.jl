@generated function struct_as_namedtuple(st)
  A = ( Expr(:(=), n, :(st.$n)) for n in fieldnames(st))
  Expr(:tuple, A...)
end

remaker_of(prob::T) where {T} = parameterless_type(T){isinplace(prob)}

# Define `remaker_of` for the types that does not (make sense to)
# implement `isinplace` trait:
for T in [
    NoiseProblem,
    SplitFunction,  # TODO: use isinplace path for type-stability
    TwoPointBVPFunction,
    ]
  @eval remaker_of(::$T) = $T
end

"""
    remake(thing; <keyword arguments>)

Re-construct `thing` with new field values specified by the keyword
arguments.
"""
function remake(thing; kwargs...)
  T = remaker_of(thing)
  T(; struct_as_namedtuple(thing)...,kwargs...)
end

isrecompile(prob::ODEProblem{iip}) where {iip} = !(typeof(prob.f) <: FunctionWrapper)

function remake(thing::ODEProblem; kwargs...)
  T = remaker_of(thing)
  tup = (struct_as_namedtuple(thing)...,kwargs...)
  if !isrecompile(thing)
    ODEProblem{isinplace(thing),false}(; tup...)
  else
    T(; struct_as_namedtuple(thing)...,kwargs...)
  end
end

function remake(thing::AbstractJumpProblem; kwargs...)
  parameterless_type(thing)(remake(thing.prob;kwargs...))
end
