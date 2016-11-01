### Basic Functionality

getindex{s}(p::ParameterizedFunction,::Val{s}) = getfield(p,s) ## Val for type-stability
