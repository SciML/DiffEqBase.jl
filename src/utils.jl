"""
`numparameters(f)`

Returns the number of parameters of `f` for the method which has the most parameters.
"""
function numargs(f)
  typ = Tuple{Any, Val{:analytic}, Vararg}
  typ2 = Tuple{Any, Type{Val{:analytic}}, Vararg} # This one is required for overloaded types
  typ3 = Tuple{Any, Val{:jac}, Vararg}
  typ4 = Tuple{Any, Type{Val{:jac}}, Vararg} # This one is required for overloaded types
  typ5 = Tuple{Any, Val{:tgrad}, Vararg}
  typ6 = Tuple{Any, Type{Val{:tgrad}}, Vararg} # This one is required for overloaded types
  numparam = maximum([(m.sig<:typ || m.sig<:typ2 || m.sig<:typ3 || m.sig<:typ4 || m.sig<:typ5 || m.sig<:typ6) ? 0 : num_types_in_tuple(m.sig) for m in methods(f)])
  return (numparam-1) #-1 in v0.5 since it adds f as the first parameter
end

function num_types_in_tuple(sig)
  length(sig.parameters)
end

function num_types_in_tuple(sig::UnionAll)
  length(Base.unwrap_unionall(sig).parameters)
end

function isinplace(f,inplace_param_number)
  numargs(f)>=inplace_param_number
end

isinplace(f::AbstractDiffEqFunction{iip}) where {iip} = iip
isinplace(f::AbstractDiffEqFunction{iip}, inplace_param_number) where {iip} = iip

macro CSI_str(str)
    return :(string("\x1b[", $(esc(str)), "m"))
end

const TYPE_COLOR = CSI"36"
const NO_COLOR = CSI"0"

macro def(name, definition)
    return quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end

using Compat.TypeUtils: typename

parameterless_type(T::Type) = typename(T).wrapper
parameterless_type(x) = parameterless_type(typeof(x))

# support functions
export check_keywords, warn_compat
function check_keywords(alg, kwargs, warnlist)
    flg = false
    for (kw, val) in kwargs
        if kw in warnlist
            if val != nothing
                flg = true
                @warn(string("The ", kw, " argument is ignored by ", alg, "."))
            end
        end
    end
    flg
end
warn_compat() =
    @warn("Please see http://docs.juliadiffeq.org/latest/basics/compatibility_chart.html")


"""
    @add_kwonly function_definition

Define keyword-only version of the `function_definition`.

    @add_kwonly function f(a, b; c=1, d=2)
        ...
    end

expands to:

    function f(x; y=1)
        ...
    end
    function f(; x = error("No argument x"), y=1)
        ...
    end
"""
macro add_kwonly(ex)
  esc(add_kwonly(ex))
end

add_kwonly(ex::Expr) = add_kwonly(Val{ex.head}, ex)

function add_kwonly(::Type{<: Val}, ex)
  error("add_only does not work with expression $(ex.head)")
end

function add_kwonly(::Union{Type{Val{:function}},
                            Type{Val{:(=)}}}, ex::Expr)
  body = ex.args[2:end]  # function body
  default_call = ex.args[1]  # e.g., :(f(a, b=2; c=3))
  kwonly_call = add_kwonly(default_call)
  if kwonly_call === nothing
    return ex
  end

  return quote
    begin
      $ex
      $(Expr(ex.head, kwonly_call, body...))
    end
  end
end

function add_kwonly(::Type{Val{:where}}, ex::Expr)
  default_call = ex.args[1]
  rest = ex.args[2:end]
  kwonly_call = add_kwonly(default_call)
  if kwonly_call === nothing
    return nothing
  end
  return Expr(:where, kwonly_call, rest...)
end

function add_kwonly(::Type{Val{:call}}, default_call::Expr)
  # default_call is, e.g., :(f(a, b=2; c=3))
  funcname = default_call.args[1]  # e.g., :f
  required = []  # required positional arguments; e.g., [:a]
  optional = []  # optional positional arguments; e.g., [:(b=2)]
  default_kwargs = []
  for arg in default_call.args[2:end]
    if isa(arg, Symbol)
      push!(required, arg)
    elseif arg.head == :(::)
      push!(required, arg)
    elseif arg.head == :kw
      push!(optional, arg)
    elseif arg.head == :parameters
      @assert default_kwargs == []  # can I have :parameters twice?
      default_kwargs = arg.args
    else
      error("Not expecting to see: $arg")
    end
  end
  if isempty(required) && isempty(optional)
    # If the function is already keyword-only, do nothing:
    return nothing
  end
  if isempty(required)
    # It's not clear what should be done.  Let's not support it at
    # the moment:
    error("At least one positional mandatory argument is required.")
  end

  kwonly_kwargs = Expr(:parameters, [
    Expr(:kw, pa, :(error($("No argument $pa"))))
    for pa in required
  ]..., optional..., default_kwargs...)
  kwonly_call = Expr(:call, funcname, kwonly_kwargs)
  # e.g., :(f(; a=error(...), b=error(...), c=1, d=2))

  return kwonly_call
end

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

"""
    undefined_exports(mod)

List symbols `export`'ed but not actually defined.
"""
function undefined_exports(mod)
  undefined = []
  for name in names(mod)
    if ! isdefined(mod, name)
      push!(undefined, name)
    end
  end
  return undefined
end
