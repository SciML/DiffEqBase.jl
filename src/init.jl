value(x) = x
promote_tspan(u0, p, tspan, prob, kwargs) = _promote_tspan(tspan, kwargs)
function _promote_tspan(tspan, kwargs)
    if (dt = get(kwargs, :dt, nothing)) !== nothing
        tspan1, tspan2, _ = promote(tspan..., dt)
        return (tspan1, tspan2)
    else
        return tspan
    end
end
isdistribution(u0) = false

function SciMLBase.tmap(args...)
    error("Zygote must be added to differentiate Zygote? If you see this error, report it.")
end

function __init__()
    @static if !isdefined(Base, :get_extension)
        @require Measurements="eff96d63-e80a-5855-80a2-b1b0885c5ab7" begin include("../ext/DiffEqBaseMeasurementsExt.jl") end

        @require MonteCarloMeasurements="0987c9cc-fe09-11e8-30f0-b96dd679fdca" begin include("../ext/DiffEqBaseMonteCarloMeasurementsExt.jl") end

        @require Unitful="1986cc42-f94f-5a68-af5c-568840ba703d" begin include("../ext/DiffEqBaseUnitfulExt.jl") end

        @require GeneralizedGenerated="6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb" begin include("../ext/DiffEqBaseGeneralizedGeneratedExt.jl") end

        @require Tracker="9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c" begin include("../ext/DiffEqBaseTrackerExt.jl") end

        @require ReverseDiff="37e2e3b7-166d-5795-8a7a-e32c996b4267" begin include("../ext/DiffEqBaseReverseDiffExt.jl") end

        @require Zygote="e88e6eb3-aa80-5325-afca-941959d7151f" begin include("../ext/DiffEqBaseZygoteExt.jl") end
        
        @require MPI="da04e1cc-30fd-572f-bb4f-1f8673147195" begin include("../ext/DiffEqBaseMPIExt.jl") end
    end
end
