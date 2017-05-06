"""
`numparameters(f)`

Returns the number of parameters of `f` for the method which has the most parameters.
"""
function numargs(f)
  typ = Tuple{Any, Val{:analytic}, Vararg}
  typ2 = Tuple{Any, Type{Val{:analytic}}, Vararg} # This one is required for overloaded types
  numparam = maximum([(m.sig<:typ || m.sig<:typ2) ? 0 : length(m.sig.parameters) for m in methods(f)]) #in v0.5, all are generic
  return (numparam-1) #-1 in v0.5 since it adds f as the first parameter
end

function isinplace(f,inplace_param_number)
  numargs(f)>=inplace_param_number
end

function isinplace{iip}(f::AbstractParameterizedFunction{iip},inplace_param_number)
  iip
end

macro def(name, definition)
    return quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end

macro muladd(ex)
  esc(to_muladd(ex))
end

function to_muladd(ex)
  is_add_operation(ex) || return ex

  all_operands = ex.args[2:end]
  mul_operands = filter(is_mul_operation, all_operands)
  odd_operands = filter(x->!is_mul_operation(x), all_operands)

  muladd_operands = collect(zip(
    to_muladd.((x->x.args[2]).(mul_operands)),
    to_muladd.((x->x.args[3]).(mul_operands))))

  if isempty(odd_operands)
    to_be_muladded = muladd_operands[1:end-1]
    last_operation = :($(muladd_operands[end][1]) * $(muladd_operands[end][2]))
  else
    to_be_muladded = muladd_operands
    last_operation = make_addition(odd_operands)
  end

  foldr(last_operation, to_be_muladded) do xs, r
    :($(Base.muladd)($(xs[1]), $(xs[2]), $r))
  end
end

is_operation(ex::Expr, op::Symbol) = ex.head == :call && !isempty(ex.args) && ex.args[1] == op
is_operation(ex, op::Symbol) = false

is_add_operation(ex) = is_operation(ex, :+)
is_mul_operation(ex) = is_operation(ex, :*)

make_addition(args) = length(args) == 1 ? args[1] : Expr(:call, :+, args...)

realtype{T}(::Type{T}) = T
realtype{T}(::Type{Complex{T}}) = T
