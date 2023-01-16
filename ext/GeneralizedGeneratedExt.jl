module GeneralizedGeneratedExt

using DiffEqBase
isdefined(Base, :get_extension) ? (using GeneralizedGenerated) :
(using ..GeneralizedGenerated)

function SciMLBase.numargs(::GeneralizedGenerated.RuntimeFn{Args}) where {Args}
    GeneralizedGenerated.from_type(Args) |> length
end

end
