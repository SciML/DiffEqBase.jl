struct OrdinaryDiffEqTag end

const NORECOMPILE_ARGUMENT_MESSAGE = """
                                     No-recompile mode is only supported for state arguments
                                     of type `Vector{Float64}`, time arguments of `Float64`
                                     and parameter arguments of type `Vector{Float64}` or
                                     `SciMLBase.NullParameters`.
                                     """

struct NoRecompileArgumentError <: Exception
    args::Any
end

function Base.showerror(io::IO, e::NoRecompileArgumentError)
    println(io, NORECOMPILE_ARGUMENT_MESSAGE)
    print(io, "Attempted arguments: ")
    print(io, e.args)
end

function unwrap_fw(fw::FunctionWrapper)
    fw.obj[]
end

# Default dispatch assumes no ForwardDiff, gets added in the new dispatch
function wrapfun_iip(ff, inputs)
    FunctionWrappersWrappers.FunctionWrappersWrapper(Void(ff), (typeof(inputs),), (Nothing,))
end

function wrapfun_oop(ff, inputs)
    FunctionWrappersWrappers.FunctionWrappersWrapper(ff, (typeof(inputs),), (typeof(inputs[1]),))
end