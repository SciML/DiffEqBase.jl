using LinearAlgebra

isscalartype(::Type{<:AbstractArray}) = false
isscalartype(::Type) = true

isfasttype(::Type{<:SubArray{<:Any,<:Any,T}}) where {T} = isfasttype(T)
isfasttype(::Type{<:Adjoint{<:Any,T}}) where {T} = isfasttype(T)
isfasttype(::Type{<:Transpose{<:Any,T}}) where {T} = isfasttype(T)
isfasttype(::Type{<:AbstractArray}) = false
isfasttype(::Type{<:Array}) = true
isfasttype(::Type) = true

@generated isfastloopable(args...) = all(isfasttype, args)

allequal(f) = true
allequal(f, x) = true
allequal(f, x, y, z...) = f(x, y) && allequal(f, y,z...)

@noinline throwdm() = throw(DimensionMismatch())

@inline @generated function fbcast(f, X, args...)
  checks = [:(axes(X))]
  if isscalartype(X)
    error("Output must be an array")
  end
  indexer = map(enumerate(args)) do ((i, arg))
    if isscalartype(arg)
      :(args[$i])
    else
      push!(checks, :(axes(args[$i])))
      :(args[$i][j])
    end
  end

  if length(checks) > 1
    loopvar = :(CartesianIndices($(first(checks))))
    quote
        allequal(===, $(checks...)) || throwdm()
        @simd ivdep for j in $loopvar
          @inbounds X[j] = f($(indexer...))
        end
        X
      end
  else
    :(fill!(X, f(args...)))
  end
end

# expression which is a call tree
# collect all the variables,
# make a function with those variables as input, and the expr as body
# call fbcast with the function and the variables list

findinputvars(vars, ex::Symbol) = push!(vars, ex=>ex)
findinputvars(vars, ex) = nothing
function findinputvars(vars, ex::Expr)
  if (ex.head in (:ref, :(.), :tuple, :macrocall))
    tmpvar = gensym("tmp")
    push!(vars, tmpvar=>ex)
  elseif ex.head == :call
    ex.args[1] == :Ref && error("Ref in fastbcast not supported, please use (a,) instead of Ref(a).")
    foreach(arg->findinputvars(vars, arg), ex.args[2:end])
  else
    error("Invalid @..")
  end
  vars
end

findinputvars(ex) = findinputvars([], ex)

macro ..(expr)
  if expr.head === :(=)
    x = expr.args[1]
    _vars = findinputvars(expr.args[2])
    _vars === nothing && @goto BC
    vars = unique(first, _vars)

    fn = quote
      @inline function ($(map(first, vars)...),)
        $(expr.args[2])
      end
    end

    quote
      if DiffEqBase.isfastloopable($x, $(map(last, vars)...))
        DiffEqBase.fbcast($fn, $x, $(map(last, vars)...))
      else
        @. $expr
      end
    end
  else
    @label BC
    :(@. $expr)
  end |> esc
end
