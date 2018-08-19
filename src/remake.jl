@generated function struct_as_namedtuple(st)
  A = ( Expr(:(=), n, :(st.$n)) for n in fieldnames(st))
  Expr(:tuple, A...)
end

"""
    remake(thing; <keyword arguments>)

Re-construct `thing` with new field values specified by the keyword
arguments.
"""
function remake(thing; kwargs...)
  T = parameterless_type(typeof(thing))
  T{isinplace(thing)}(; struct_as_namedtuple(thing)...,kwargs...)
end
