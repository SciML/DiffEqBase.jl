# the following are setup per how integrators are implemented in OrdinaryDiffEq and
# StochasticDiffEq and provide dispatch points that JumpProcesses and others can use.

get_tstops(integ) = integ.opts.tstops
get_tstops_array(integ) = get_tstops(integ).valtree
get_max_tstops(integ) = maximum(get_tstops_array(integ))
