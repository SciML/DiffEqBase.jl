struct StandardBVProblem end

struct BVProblem{uType,tType,isinplace,P,F,bF,PT,CB} <: AbstractBVProblem{uType,tType,isinplace}
    f::F
    bc::bF
    u0::uType
    tspan::tType
    p::P
    problem_type::PT
    callback::CB
    @add_kwonly function BVProblem{iip}(f::AbstractODEFunction,bc,u0,tspan,p=NullParameters(),
                            problem_type=StandardBVProblem();
                            callback=nothing) where {iip}
        _tspan = promote_tspan(tspan)
        new{typeof(u0),typeof(tspan),iip,typeof(p),
                  typeof(f),typeof(bc),
                  typeof(problem_type),typeof(callback)}(
                  f,bc,u0,_tspan,p,
                  problem_type,callback)
    end

    function BVProblem{iip}(f,bc,u0,tspan,p=NullParameters();kwargs...) where {iip}
        BVProblem(convert(ODEFunction{iip},f),bc,u0,tspan,p;kwargs...)
      end
end

function BVProblem(f::AbstractODEFunction,bc,u0,tspan,args...;kwargs...)
    BVProblem{DiffEqBase.isinplace(f,4)}(f,bc,u0,tspan,args...;kwargs...)
end

function BVProblem(f,bc,u0,tspan,p=NullParameters();kwargs...)
    BVProblem(convert(ODEFunction,f),bc,u0,tspan,p;kwargs...)
end

# convenience interfaces:
# Allow any previous timeseries solution
function BVProblem(f,bc,sol::T,tspan,p=NullParameters();kwargs...) where {T<:AbstractTimeseriesSolution}
    BVProblem(f,bc,sol.u,tspan,p)
end
# Allow a function of time for the initial guess
function BVProblem(f,bc,initialGuess::T,tspan::AbstractVector,p=NullParameters();kwargs...) where {T}
    u0 = [ initialGuess( i ) for i in tspan]
    BVProblem(f,bc,u0,(tspan[1],tspan[end]),p)
end

struct TwoPointBVPFunction{bF}
    bc::bF
end
TwoPointBVPFunction(; bc = error("No argument bc")) = TwoPointBVPFunction(bc)
(f::TwoPointBVPFunction)(residual, ua, ub, p) = f.bc(residual, ua, ub, p)
(f::TwoPointBVPFunction)(residual, u, p) = f.bc(residual, u[1], u[end], p)

struct TwoPointBVProblem{iip} end
function TwoPointBVProblem(f,bc,u0,tspan,p=NullParameters();kwargs...)
    iip = DiffEqBase.isinplace(f,4)
    TwoPointBVProblem{iip}(f,bc,u0,tspan,p;kwargs...)
end
function TwoPointBVProblem{iip}(f,bc,u0,tspan,p=NullParameters();kwargs...) where {iip}
    BVProblem{iip}(f,TwoPointBVPFunction(bc),u0,tspan,p;kwargs...)
end
