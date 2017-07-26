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

"""
    @muladd ex

Convert every combined multiplication and addition in `ex` into a call of `muladd`. If any
of the involved operators or operands is dotted, `muladd` is applied as a "dot call".
"""
macro muladd(ex)
  esc(to_muladd(ex))
end

function to_muladd(ex::Expr)
    if !is_add_operation(ex)
        if ex.head == :macrocall
            # expand macros first (enables use of @. inside of @muladd expression)
            return to_muladd(macroexpand(ex))
        else
            # if expression is no sum apply the reduction to its arguments
            return Expr(ex.head, to_muladd.(ex.args)...)
        end
    end

    # retrieve summands of addition and split them into two groups, one with expressions
    # of multiplications and one with other expressions
    all_operands = to_muladd.(operands(ex))
    mul_operands = filter(is_mul_operation, all_operands)
    odd_operands = filter(x->!is_mul_operation(x), all_operands)

    # define summands that are reduced with muladd and the initial element of the reduction
    if isempty(odd_operands)
        # if all summands are multiplications one of these summands is
        # the initial element of the reduction
        to_be_muladded = mul_operands[1:end-1]
        last_operation = mul_operands[end]
    else
        to_be_muladded = mul_operands

        # expressions that are no multiplications are summed up in a separate expression
        # that is the initial element of the reduction
        # if the original addition was a dot call this expression also is a dot call
        if length(odd_operands) == 1
            last_operation = odd_operands[1]
        elseif isdotcall(ex)
            last_operation = Expr(:., :+, Expr(:tuple, odd_operands...))
        else
            last_operation = Expr(:call, :+, odd_operands...)
        end
    end

    # reduce sum to a composition of muladd
    foldr(last_operation, to_be_muladded) do xs, r
        # retrieve factors of multiplication that will be reduced next
        xs_operands = operands(xs)

        # first factor is always first operand
        xs_factor1 = xs_operands[1]

        # second factor is an expression of a multiplication if there are more than
        # two operands
        # if the original multiplication was a dot call this expression also is a dot call
        if length(xs_operands) == 2
            xs_factor2 = xs_operands[2]
        elseif isdotcall(xs)
            xs_factor2 = Expr(:., :*, Expr(:tuple, xs_operands[2:end]...))
        else
            xs_factor2 = Expr(:call, :*, xs_operands[2:end]...)
        end

        # create a dot call if any of the involved operators or operands is a dot call
        if any(isdotcall, (ex, xs, xs_factor1, xs_factor2, r))
            Expr(:., Base.muladd, Expr(:tuple, xs_factor1, xs_factor2, r))
        else
            Expr(:call, Base.muladd, xs_factor1, xs_factor2, r)
        end
    end
end
to_muladd(ex) = ex

"""
    isoperation(ex, op::Symbol)

Determine whether `ex` is a call or "dot call" of operation `op`.
"""
isoperation(ex::Expr, op::Symbol) = !isempty(ex.args) &&
    ((ex.head == :call && (ex.args[1] == op || ex.args[1] == Symbol('.', op))) ||
     (ex.head == :. && ex.args[1] == op))
isoperation(ex, op::Symbol) = false

is_add_operation(ex) = isoperation(ex, :+)
is_mul_operation(ex) = isoperation(ex, :*)

"""
    isdotcall(ex)

Determine whether `ex` is a dot call.
"""
isdotcall(ex::Expr) = ex.head == :. || (ex.head == :call && !isempty(ex.args) &&
                                        first(string(ex.args[1])) == '.')
isdotcall(ex) = false

"""
    operands(ex)

Return arguments of function call in `ex`.
"""
function operands(ex::Expr)
    if ex.head == :. && length(ex.args) == 2 && typeof(ex.args[2]) <: Expr
        ex.args[2].args
    else
        ex.args[2:end]
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
