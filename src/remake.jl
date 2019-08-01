@generated function struct_as_namedtuple(st)
  A = ( Expr(:(=), n, :(st.$n)) for n in fieldnames(st))
  Expr(:tuple, A...)
end

Base.@pure remaker_of(prob::T) where {T} = parameterless_type(T){isinplace(prob)}

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

function remake(thing::AbstractDiffEqFunction{iip}; kwargs...) where iip
  T = ODEFunction{iip}
  T(; struct_as_namedtuple(thing)...,kwargs...)
end

isrecompile(prob::ODEProblem{iip}) where {iip} = (prob.f isa ODEFunction) ? !(typeof(prob.f.f) <: FunctionWrapper) : true

function remake(thing::ODEProblem; kwargs...)
  T = remaker_of(thing)
  tup = merge(struct_as_namedtuple(thing),kwargs)
  if !isrecompile(thing)
    if isinplace(thing)
      f = wrapfun_iip(unwrap_fw(tup.f.f),(tup.u0,tup.u0,tup.p,tup.tspan[1]))
    else
      f = wrapfun_oop(unwrap_fw(tup.f.f),(tup.u0,tup.p,tup.tspan[1]))
    end
    tup2 = (f = convert(ODEFunction{isinplace(thing)},f),)
    tup = merge(tup, tup2)
  end
  T(; tup...)
end

function remake(thing::AbstractJumpProblem; kwargs...)
  parameterless_type(thing)(remake(thing.prob;kwargs...))
end
