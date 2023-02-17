module DiffEqBaseDistributionsExt

using Distributions, DiffEqBase

DiffEqBase.handle_distribution_u0(_u0::Distributions.Sampleable) = rand(_u0)
DiffEqBase.isdistribution(_u0::Distributions.Sampleable) = true

end
