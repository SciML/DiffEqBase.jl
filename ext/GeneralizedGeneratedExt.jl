module GeneralizedGeneratedExt

using GeneralizedGenerated, DiffEqBase

function SciMLBase.numargs(::GeneralizedGenerated.RuntimeFn{Args}) where {Args}
    GeneralizedGenerated.from_type(Args) |> length
end

end