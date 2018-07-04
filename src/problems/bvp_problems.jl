struct StandardBVProblem end

struct BVProblem{uType,tType,isinplace,P,J,F,bF,PT,CB,MM} <: AbstractBVProblem{uType,tType,isinplace}
    f::F
    bc::bF
    u0::uType
    tspan::tType
    p::P
    jac_prototype::J
    problem_type::PT
    callback::CB
    mass_matrix::MM
    @add_kwonly function BVProblem{iip}(f,bc,u0,tspan,p=nothing,
                            problem_type=StandardBVProblem();
                            jac_prototype = nothing,
                            callback=nothing,mass_matrix=I) where {iip}
        _tspan = promote_tspan(tspan)
        new{typeof(u0),typeof(tspan),iip,typeof(p),typeof(jac_prototype),
                  typeof(f),typeof(bc),
                  typeof(problem_type),typeof(callback),typeof(mass_matrix)}(
                  f,bc,u0,_tspan,p,jac_prototype,
                  problem_type,callback,mass_matrix)
    end
end

function BVProblem(f,bc,u0::AbstractArray,tspan,p=nothing;kwargs...)
    iip = DiffEqBase.isinplace(f,4)
    BVProblem{iip}(f,bc,u0,tspan,p;kwargs...)
end

# convenience interfaces:
# Allow any previous timeseries solution
function BVProblem(f,bc,sol::T,tspan,p=nothing;kwargs...) where {T<:AbstractTimeseriesSolution}
    BVProblem(f,bc,sol.u,tspan,p)
end
# Allow a function of time for the initial guess
function BVProblem(f,bc,initialGuess::T,tspan::AbstractVector,p=nothing;kwargs...) where {T}
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
function TwoPointBVProblem(f,bc,u0,tspan,p=nothing;kwargs...)
    iip = DiffEqBase.isinplace(f,4)
    TwoPointBVProblem{iip}(f,bc,u0,tspan,p;kwargs...)
end
function TwoPointBVProblem{iip}(f,bc,u0,tspan,p=nothing;kwargs...) where {iip}
    BVProblem{iip}(f,TwoPointBVPFunction(bc),u0,tspan,p;kwargs...)
end
