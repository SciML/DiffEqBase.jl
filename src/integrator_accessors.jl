# the following are setup per how integrators are implemented in OrdinaryDiffEq and
# StochasticDiffEq and provide dispatch points that JumpProcesses and others can use.

get_tstops(integ::DEIntegrator) =
    error("get_tstops not implemented for integrators of type $(typeof(integ))")
get_tstops_array(integ::DEIntegrator) =
    error("get_tstops_array not implemented for integrators of type $(typeof(integ))")
get_tstops_max(integ::DEIntegrator) =
    error("get_tstops_max not implemented for integrators of type $(typeof(integ))")
