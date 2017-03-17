__precompile__()
module TypeMacros
export @abstract_type, @struct_type, @public_abstract_type, @public_abstract_subs

macro abstract_type(exp) ; _at(exp) ; end
macro struct_type(exp)   ; _st(exp) ; end
macro public_abstract_type(exp); _pat(exp); end
macro public_abstract_subs(exp); _pas(exp); end

@static if VERSION < v"0.6.0-dev.2746"
    include("typemacros5.jl")
else
    include("typemacros6.jl")
end

_type_def_error() = throw(ArgumentError("not a valid type definition"))
function _at(exp)
    isa(exp, Symbol) && return _qa(exp)
    isa(exp, Expr) || _type_def_error()
    if exp.head == :curly
        return _qa(exp)
    elseif exp.head == :(<:)
        return _qa(exp.args[1], exp.args[2])
    else
        _type_def_error()
    end
end
    
function _st(exp)
    isa(exp, Symbol) && return _qs(exp)
    isa(exp, Expr) || _type_def_error()
    if exp.head == :curly
        return _qs(exp)
    elseif exp.head == :(<:)
        return _qs(exp.args[1], exp.args[2])
    else
        _type_def_error()
    end
end

function _pat(exp)
    blk = quote end
    if isa(exp, Symbol)
        push!(blk.args, _qa(exp))
        sym = exp
    elseif !isa(exp, Expr)
        _type_def_error()
    elseif exp.head == :curly
        push!(blk.args, _qa(exp))
        sym = exp.args[1]
    elseif exp.head == :(<:)
        push!(blk.args, _qa(exp.args[1], exp.args[2]))
        exp1 = exp.args[1]
        if isa(exp1, Symbol)
            sym = exp1
        elseif isa(exp1, Expr) && exp1.head == :curly
            sym = exp1.args[1]
        else
            _type_def_error()
        end
    else
        _type_def_error()
    end
    push!(blk.args, :( export $(esc(sym)) ))
    push!(blk.args, :nothing)
    blk.head = :toplevel
    blk
end

function _pas(exp)
    (isa(exp, Expr) && exp.head == :(<:)) || _type_def_error()
    length(exp.args) == 2 || _type_def_error()
    (isa(exp.args[1], Expr) && exp.args[1].head == :tuple) || _type_def_error()
    blk = quote end
    syms = exp.args[1].args
    sup  = exp.args[2]
    for s in syms
        if isa(s, Symbol)
            push!(blk.args, _qa(s, sup))
            push!(blk.args, :( export $(esc(s)) ))
        elseif isa(s, Expr) && s.head == :curly
            push!(blk.args, _qa(s, sup))
            push!(blk.args, :( export $(esc(s.args[1])) ) )
        else
            _type_def_error()
        end
    end
    push!(blk.args, :nothing)
    blk.head = :toplevel
    blk
end
end
