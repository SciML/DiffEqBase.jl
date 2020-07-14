"""
$(TYPEDEF)

Holds a tableau which defines an explicit Runge-Kutta method.
"""
struct RKTableau{AType,bType,bbType,cType,
                 fsal,explicit,stages,order,adaptiveorder,name} <: ODERKTableau
  A::AType
  c::cType
  b::bType
  bEEst::bbType
end
function RKTableau(A, c, b, order; adaptiveorder=0, bEEst=nothing, name=:RK, kwargs...)
  stages = length(b)
  eeststages = bEEst === nothing ? stages : length(bEEst)
  size(A, 1) == size(A, 2) == stages == eeststages || error("The size of A is not consistent with b or embedded b.")
  isfsal = iszero(c[1]) && isone(c[stages]) && A[stages, :] == b
  isexplicit = true
  for j in axes(A, 2), i in j:stages
    isexplicit &= iszero(A[i, j])
  end

  RKTableau{typeof(A), typeof(b), typeof(bEEst), typeof(c),
            isfsal, isexplicit, stages, order, adaptiveorder, name}(
    A, c, b, bEEst
  )
end

const ExplicitRKTableau{AType,bType,bbType,cType,fsal,stages,order,adaptiveorder} = RKTableau{AType,bType,bbType,cType,fsal,true,stages,order,adaptiveorder}
const ImplicitRKTableau{AType,bType,bbType,cType,fsal,stages,order,adaptiveorder} = RKTableau{AType,bType,bbType,cType,fsal,false,stages,order,adaptiveorder}

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

for f in [:isexplicit, :alg_order, :alg_adaptive_order, :isadaptive, :isfsal, :alg_stages, :get_tableau_A, :get_tableau_b, :get_tableau_bEEst, :get_tableau_c]
  @eval $f(::T) where {T<:Any} = error("`$f(::$T)` is not defined.")
end
