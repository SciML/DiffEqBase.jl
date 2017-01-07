"""
`numparameters(f)`

Returns the number of parameters of `f` for the method which has the most parameters.
"""
function numparameters(f)
  numparam = maximum([length(m.sig.parameters) for m in methods(f)]) #in v0.5, all are generic
  return (numparam-1) #-1 in v0.5 since it adds f as the first parameter
end

function isinplace(f,inplace_param_number)
  numparameters(f)>=inplace_param_number
end

function isinplace{iip}(f::AbstractParameterizedFunction{iip},inplace_param_number)
  iip
end

macro def(name, definition)
    return quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end
