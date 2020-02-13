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

    push!(value_types, timestamp_type)
    push!(value_names, :timestamp)

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

    si = DESolutionIterator{NamedTuple{(value_names...,),Tuple{value_types...}},typeof(sol)}(sol)

    return si
end

function Base.length(iter::DESolutionIterator)
    return length(iter.sol)
end

Base.eltype(::Type{DESolutionIterator{T,TS}}) where {T,TS} = T

@generated function Base.iterate(iter::DESolutionIterator{T,S}, state=1) where {T,S}
    columns = []
    push!(columns, :(iter.sol.t[state]))
    nfield = length(fieldnames(T))
    for i in 1:nfield-1
        if nfield > 2
            push!(columns, :(iter.sol.u[state][$i]))
        else
            push!(columns, :(iter.sol.u[state]))
        end
    end

    quote
        if state > length(iter)
            return nothing
        else
            return $(T)(($(columns...),)), state+1
        end
    end
end
