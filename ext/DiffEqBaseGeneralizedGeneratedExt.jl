module DiffEqBaseGeneralizedGeneratedExt

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
