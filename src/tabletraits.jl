struct DESolutionIterator{T,S}
    sol::S
end

IteratorInterfaceExtensions.isiterable(sol::DESolution) = true
TableTraits.isiterabletable(sol::DESolution) = true

function IteratorInterfaceExtensions.getiterator(sol::DESolution)
    timestamp_type = eltype(sol.t)
    value_type = eltype(sol.u)

    value_types = Type[]
    value_names = Symbol[]

    if value_type<:AbstractArray
        for i in 1:length(sol.u[1])
            push!(value_types, eltype(sol.u[1]))
            if has_syms(sol.prob.f)
                push!(value_names, sol.prob.f.syms[i])
            else
                push!(value_names, Symbol("value$i"))
            end
          end
    else
        push!(value_types, value_type)
        if has_syms(sol.prob.f)
            push!(value_names, sol.prob.f.syms[1])
        else
            push!(value_names, :value)
        end
    end

    col_expressions = Array{Expr,1}()

    # Add timestamp column
    push!(col_expressions, Expr(:(::), :timestamp, timestamp_type))

    for i in 1:length(value_types)
        push!(col_expressions, Expr(:(::), value_names[i], value_types[i]))
    end

    t_expr = NamedTuples.make_tuple(col_expressions)

    t2 = :(DESolutionIterator{Float64,Float64})
    t2.args[2] = t_expr
    t2.args[3] = typeof(sol)

    t = eval(t2)

    si = t(sol)

    return si
end

function Base.length(iter::DESolutionIterator)
    return length(iter.sol)
end

function Base.eltype(iter::DESolutionIterator{T,S}) where {T,S}
    return T
end

Base.eltype(::Type{DESolutionIterator{T,TS}}) where {T,TS} = T

function Base.start(iter::DESolutionIterator)
    return 1
end

@generated function Base.next(iter::DESolutionIterator{T,S}, state) where {T,S}
    constructor_call = Expr(:call, :($T))
    push!(constructor_call.args, :(iter.sol.t[i]))
    for i in 1:length(T.parameters)-1
        if length(T.parameters)>2
            push!(constructor_call.args, :(iter.sol.u[i][$i]))
        else
            push!(constructor_call.args, :(iter.sol.u[i]))
        end
    end

    quote
        i = state
        a = $constructor_call
        return a, state+1
    end
end

function Base.done(iter::DESolutionIterator, state)
    return state>length(iter.sol)
end
