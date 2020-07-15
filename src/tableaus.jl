"""
  $(TYPEDEF)

Holds a tableau which defines a Runge-Kutta method.
"""
struct RKTableau{AType,bType,bbType,cType,
                 fsal,explicit,stages,order,adaptiveorder,name} <: ODERKTableau
  A::AType
  c::cType
  b::bType
  bEEst::bbType
end
function RKTableau(A, c, b, order; adaptiveorder=0, αEEst=nothing, bEEst=αEEst, name=:RK, kwargs...)
  #TODO
  αEEst !== nothing && Base.depwarn("RKTableau(A, c, b, order; αEEst=X) is deprecated. Please use RKTableau(A, c, b, order; bEEst=X) instead.", ((Base.Core).Typeof(RKTableau)).name.mt.name)
  stages = length(b)
  bEEst = bEEst === nothing ? b : bEEst
  eeststages = length(bEEst)
  size(A, 1) == size(A, 2) == stages == eeststages || error("The size of A is not consistent with b or embedded b.")
  isfsal = iszero(c[1]) && isone(c[stages]) && A[stages, :] == b
  isexplicit = istril(A, -1) && iszero(c[1])

  RKTableau{typeof(A), typeof(b), typeof(bEEst), typeof(c),
            isfsal, isexplicit, stages, order, adaptiveorder, name}(
    A, c, b, bEEst
  )
end

const ExplicitRKTableau{AType,bType,bbType,cType,fsal,stages,order,adaptiveorder} = RKTableau{AType,bType,bbType,cType,fsal,true,stages,order,adaptiveorder}
const ImplicitRKTableau{AType,bType,bbType,cType,fsal,stages,order,adaptiveorder} = RKTableau{AType,bType,bbType,cType,fsal,false,stages,order,adaptiveorder}

###
### Deprecation
###

Base.@deprecate ExplicitRKTableau(args...; kw...) RKTableau(args...; kw...) false
Base.@deprecate ImplicitRKTableau(args...; kw...) RKTableau(args...; kw...) false

const DEPWARN_COUNT = Ref(0)

function Base.getproperty(tab::RKTableau, s::Symbol)
  verbose = DEPWARN_COUNT[] <= 10
  if s === :α
    DEPWARN_COUNT[] += 1
    verbose && @warn "tab.α is deprecated. Please use tab.b instead."
    getfield(tab, :b)
  elseif s === :αEEst
    DEPWARN_COUNT[] += 1
    verbose && @warn "tab.αEEst is deprecated. Please use tab.bEEst instead."
    getfield(tab, :bEEst)
  elseif s === :order
    DEPWARN_COUNT[] += 1
    verbose && @warn "tab.order is deprecated. Please use DiffEqBase.alg_order(tab) instead."
    alg_order(tab)
  elseif s === :adaptiveorder
    DEPWARN_COUNT[] += 1
    verbose && @warn "tab.order is deprecated. Please use DiffEqBase.alg_adaptive_order(tab) instead."
    alg_adaptive_order(tab)
  elseif s === :stages
    DEPWARN_COUNT[] += 1
    verbose && @warn "tab.stage is deprecated. Please use DiffEqBase.alg_stages(tab) instead."
    alg_stages(tab)
  else
    getfield(tab, s)
  end
end

###
### Interface
###

isexplicit(tab::ExplicitRKTableau) = true
isexplicit(tab::ImplicitRKTableau) = false

alg_order(tab::RKTableau{AType,bType,bbType,cType,fsal,explicit,stages,order,adaptiveorder}) where {AType,bType,bbType,cType,fsal,explicit,stages,order,adaptiveorder} = order
alg_adaptive_order(tab::RKTableau{AType,bType,bbType,cType,fsal,explicit,stages,order,adaptiveorder}) where {AType,bType,bbType,cType,fsal,explicit,stages,order,adaptiveorder} = adaptiveorder

isadaptive(tab::ODERKTableau) = alg_adaptive_order(tab) !== 0

isfsal(tab::RKTableau{AType,bType,bbType,cType,fsal,explicit,stages,order,adaptiveorder}) where {AType,bType,bbType,cType,fsal,explicit,stages,order,adaptiveorder} = fsal

alg_stages(tab::RKTableau{AType,bType,bbType,cType,fsal,explicit,stages,order,adaptiveorder}) where {AType,bType,bbType,cType,fsal,explicit,stages,order,adaptiveorder} = stages

get_tableau_A(tab::ODERKTableau) = tab.A
get_tableau_b(tab::ODERKTableau) = tab.b
get_tableau_bEEst(tab::ODERKTableau) = tab.bEEst
get_tableau_c(tab::ODERKTableau) = tab.c
get_tableau_name(tab::RKTableau{AType,bType,bbType,cType,fsal,explicit,stages,order,adaptiveorder,name}) where {AType,bType,bbType,cType,fsal,explicit,stages,order,adaptiveorder,name} = name

for f in [:isexplicit, :alg_order, :alg_adaptive_order, :isadaptive, :isfsal, :alg_stages, :get_tableau_A, :get_tableau_b, :get_tableau_bEEst, :get_tableau_c, :get_tableau_name]
  @eval $f(::T) where {T<:Any} = error("`$f(::$T)` is not defined.")
end

###
### Pretty printing
###

function algstats(tab::ODERKTableau)
  stats = String[]
  push!(stats, isexplicit(tab) ? "explicit" : "implicit")
  isfsal(tab) && push!(stats, "FSAL")
  stages = alg_stages(tab)
  push!(stats, "stages=$stages")
  order = alg_order(tab)
  push!(stats, "order=$order")
  if isadaptive(tab)
    adaptiveorder = alg_adaptive_order(tab)
    push!(stats, "embedded order=$adaptiveorder")
  end
  return join(stats, ", ")
end

struct ExplicitA{T} <: AbstractMatrix{T}
  data::AbstractMatrix{T}
end
Base.getindex(A::ExplicitA, args...) = getindex(A.data, args...)
Base.size(A::ExplicitA, args...) = size(A.data, args...)

Base.replace_in_print_matrix(A::ExplicitA, i::Integer, j::Integer, s::AbstractString) = i > j ? s : Base.replace_with_centered_mark(s)

function Base.show(io::IO, tab::ODERKTableau)
  println(io, get_tableau_name(tab),
          "{",
          algstats(tab),
          "}",
          ":\n", "A:")
  A = get_tableau_A(tab)
  Base.print_matrix(io, isexplicit(tab) ? ExplicitA(A) : A)
  println(io, "\nb:")
  Base.print_matrix(io, get_tableau_b(tab)')
  if isadaptive(tab)
    println(io, "\nembedded b:")
    Base.print_matrix(io, get_tableau_bEEst(tab)')
  end
  println(io, "\nc:")
  Base.print_matrix(io, get_tableau_c(tab))
end
