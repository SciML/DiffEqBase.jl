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
    if !isdefined(Base, :get_extension)
        @require Measurements="eff96d63-e80a-5855-80a2-b1b0885c5ab7" begin
            include("../ext/MeasurementsExt.jl")
        end

        @require MonteCarloMeasurements="0987c9cc-fe09-11e8-30f0-b96dd679fdca" begin
            include("../ext/MonteCarloMeasurementsExt.jl")
        end

        @require Unitful="1986cc42-f94f-5a68-af5c-568840ba703d" begin
            include("../ext/UnitfulExt.jl")
        end

        @require GeneralizedGenerated="6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb" begin
            include("../ext/GeneralizedGeneratedExt.jl")
        end
    end
end
