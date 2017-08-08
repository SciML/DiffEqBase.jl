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

function isinplace(f::AbstractParameterizedFunction{iip},inplace_param_number) where iip
  iip
end

macro def(name, definition)
    return quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end

realtype{T}(::Type{T}) = T
realtype{T}(::Type{Complex{T}}) = T

using Compat.TypeUtils: typename

if :wrapper in fieldnames(TypeName)
    parameterless_type(T::Type) = typename(T).wrapper
else
    parameterless_type(T::Type) = typename(T).primary
end

parameterless_type(x) = parameterless_type(typeof(x))

# support functions
export check_keywords, warn_compat
function check_keywords(alg, kwargs, warnlist)
    flg = false
    for (kw, val) in kwargs
        if kw in warnlist
            if val != nothing
                flg = true
                warn(string("The ", kw, " argument is ignored by ", alg, "."))
            end
        end
    end
    flg
end
warn_compat() =
    warn("Please see http://docs.juliadiffeq.org/latest/basics/compatibility_chart.html")
