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
  numparm = maximum([length(m.sig.parameters) for m in methods(f)]) #in v0.5, all are generic
  return (numparm-1) #-1 in v0.5 since it adds f as the first parameter
end
