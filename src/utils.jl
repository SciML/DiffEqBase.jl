"""
`numparameters(f)`

Returns the number of parameters of `f` for the method which has the most parameters.
"""
function numparameters(f)
  numparam = maximum([length(m.sig.parameters) for m in methods(f)]) #in v0.5, all are generic
  return (numparam-1) #-1 in v0.5 since it adds f as the first parameter
end

# Check for the function overloads

jac_exists(f) = !isempty(methods(f,Tuple{Type{Val{:jac}}, Any,Any,Any}))
invjac_exists(f) = !isempty(methods(f,Tuple{Type{Val{:invjac}}, Any,Any,Any}))
hes_exists(f) = !isempty(methods(f,Tuple{Type{Val{:hes}}, Any,Any,Any}))
invhes_exists(f) = !isempty(methods(f,Tuple{Type{Val{:invhes}}, Any,Any,Any}))
paramjac_exists(f) = !isempty(methods(f,Tuple{Type{Val{:paramjac}}, Any,Any,Any,Any}))
pfunc_exists(f,p::Symbol) = !isempty(methods(f,Tuple{Type{Val{p}}, Any,Any,Any}))
pderiv_exists(f,p::Symbol) = !isempty(methods(f,Tuple{Type{Val{p}},Val{:deriv},Any,Any,Any}))
pfunc_exists(f) = !isempty(methods(f,Tuple{Type{Val{f.params[1]}}, Any,Any,Any,Any}))
pderiv_exists(f) = !isempty(methods(f,Tuple{Type{Val{f.params[1]}},Val{:deriv},Any,Any,Any}))
