Tables.isrowtable(::Type{<:DESolution}) = true
Tables.columns(x::DESolution) = Tables.columntable(Tables.rows(x))

struct DESolutionRows{T, U}
    names::Vector{Symbol}
    types::Vector{Type}
    lookup::Dict{Symbol, Int}
    t::T
    u::U
end

DESolutionRows(names, types, t, u) = DESolutionRows(names, types, Dict(nm => i for (i, nm) in enumerate(names)), t, u)

Base.length(x::DESolutionRows) = length(x.u)
Base.eltype(x::DESolutionRows{T, U}) where {T, U} = DESolutionRow{eltype(T), eltype(U)}
Base.iterate(x::DESolutionRows, st=1) = st > length(x) ? nothing : (DESolutionRow(x.names, x.lookup, x.t[st], x.u[st]), st + 1)

function Tables.rows(sol::DESolution)
    VT = eltype(sol.u)
    if VT <: AbstractArray
        N = length(sol.u[1])
        names = [:timestamp, (has_syms(sol.prob.f) ? (sol.prob.f.syms[i] for i = 1:N) : (Symbol("value", i) for i = 1:N))...]
        types = Type[eltype(sol.t), (eltype(sol.u[1]) for i = 1:N)...]
    else
        names = [:timestamp, has_syms(sol.prob.f) ? sol.prob.f.syms[1] : :value]
        types = Type[eltype(sol.t), VT]
    end
    return DESolutionRows(names, types, sol.t, sol.u)
end

Tables.schema(x::DESolutionRows) = Tables.Schema(x.names, x.types)

struct DESolutionRow{T, U} <: Tables.AbstractRow
    names::Vector{Symbol}
    lookup::Dict{Symbol, Int}
    t::T
    u::U
end

Tables.columnnames(x::DESolutionRow) = getfield(x, :names)
Tables.getcolumn(x::DESolutionRow, i::Int) = i == 1 ? getfield(x, :t) : getfield(x, :u)[i - 1]
Tables.getcolumn(x::DESolutionRow, nm::Symbol) = nm === :timestamp ? getfield(x, :t) : getfield(x, :u)[getfield(x, :lookup)[nm] - 1]

IteratorInterfaceExtensions.isiterable(sol::DESolution) = true
TableTraits.isiterabletable(sol::DESolution) = true

IteratorInterfaceExtensions.getiterator(sol::DESolution) =
    Tables.datavaluerows(Tables.rows(sol))
