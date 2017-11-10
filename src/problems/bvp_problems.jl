struct StandardBVProblem end

struct BVProblem{uType,tType,isinplace,F,bF,PT,CB,MM} <: AbstractBVProblem{uType,tType,isinplace}
    f::F
    bc::bF
    u0::uType
    tspan::Tuple{tType,tType}
    problem_type::PT
    callback::CB
    mass_matrix::MM
    function BVProblem{iip}(f,bc,u0,tspan,problem_type=StandardBVProblem();
                      callback=nothing,mass_matrix=I) where {iip}
        new{typeof(u0),eltype(tspan),iip,typeof(f),typeof(bc),
                  typeof(problem_type),typeof(callback),typeof(mass_matrix)}(
                  f,bc,u0,tspan,problem_type,callback,mass_matrix)
    end
end

function BVProblem(f,bc,u0::Array,tspan;kwargs...)
    iip = DiffEqBase.isinplace(f,3)
    BVProblem{iip}(f,bc,u0,tspan;kwargs...)
end

function BVProblem(f,bc,initialGuess::Any,tspan::StepRangeLen;kwargs...)
    u0 = [ initialGuess( i ) for i in tspan]
    BVProblem(f,bc,u0,[tspan[1],tspan[end]])
end

struct TwoPointBVPFunction{bF}
    bc::bF
end
(f::TwoPointBVPFunction)(residual, ua, ub) = f.bc(residual, ua, ub)
(f::TwoPointBVPFunction)(residual, u) = f.bc(residual, u[1], u[end])

struct TwoPointBVProblem{iip} end
function TwoPointBVProblem(f,bc,u0,tspan;kwargs...)
    iip = DiffEqBase.isinplace(f,3)
    TwoPointBVProblem{iip}(f,bc,u0,tspan;kwargs...)
end
function TwoPointBVProblem{iip}(f,bc,u0,tspan;kwargs...) where {iip}
    BVProblem{iip}(f,TwoPointBVPFunction(bc),u0,tspan;kwargs...)
end
