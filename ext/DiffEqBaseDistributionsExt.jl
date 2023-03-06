module DiffEqBaseDistributionsExt

if isdefined(Base, :get_extension)
    using Distributions
    using DiffEqBase
else
    using ..Distributions
    using ..DiffEqBase
end

DiffEqBase.handle_distribution_u0(_u0::Distributions.Sampleable) = rand(_u0)
DiffEqBase.isdistribution(_u0::Distributions.Sampleable) = true

end
