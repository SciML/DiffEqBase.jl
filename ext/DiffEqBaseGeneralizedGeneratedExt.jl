module DiffEqBaseGeneralizedGeneratedExt

if isdefined(Base, :get_extension)
    using DiffEqBase
    using GeneralizedGenerated
else
    using ..DiffEqBase
    using ..GeneralizedGenerated
end

function SciMLBase.numargs(::GeneralizedGenerated.RuntimeFn{Args}) where {Args}
    GeneralizedGenerated.from_type(Args) |> length
end

end
