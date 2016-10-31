macro def(name, definition)
    quote
        macro $name()
            esc($(Expr(:quote, definition)))
        end
    end
end

"""
`numparameters(f)`

Returns the number of parameters of `f` for the method which has the most parameters.
"""
function numparameters(f)
  #=
  if length(methods(f))>1
    warn("Number of methods for f is greater than 1. Choosing linearity based off of method with most parameters")
  end
  =#
  numparm = maximum([length(m.sig.parameters) for m in methods(f)]) #in v0.5, all are generic
  return (numparm-1) #-1 in v0.5 since it adds f as the first parameter
end
