struct OrdinaryDiffEqTag end

const dualT = ForwardDiff.Dual{ForwardDiff.Tag{OrdinaryDiffEqTag, Float64}, Float64, 1}
const NORECOMPILE_SUPPORTED_ARGS = (Tuple{Vector{Float64}, Vector{Float64},
                                          Vector{Float64}, Float64},
                                    Tuple{Vector{Float64}, Vector{Float64},
                                          SciMLBase.NullParameters, Float64})
const arglists = (Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64},
                  Tuple{Vector{Float64}, Vector{Float64}, SciMLBase.NullParameters, Float64
                        },
                  Tuple{Vector{dualT}, Vector{Float64}, Vector{Float64}, dualT},
                  Tuple{Vector{dualT}, Vector{dualT}, Vector{Float64}, Float64},
                  Tuple{Vector{dualT}, Vector{dualT}, SciMLBase.NullParameters, Float64},
                  Tuple{Vector{dualT}, Vector{Float64}, SciMLBase.NullParameters, dualT})
const iip_returnlists = ntuple(x -> Nothing, length(arglists))
function void(@nospecialize(f::Function))
    Base.Experimental.@opaque (args...)->f(args...)
end

const oop_returnlists = (Vector{Float64},Vector{Float64},
                         ntuple(x -> Vector{dualT}, length(arglists)-2)...)

function typestablemapping(@nospecialize(f::Function))
    Base.Experimental.@opaque (args...)->f(args...)
end

const NORECOMPILE_ARGUMENT_MESSAGE = """
                                     No-recompile mode is only supported for state arguments
                                     of type `Vector{Float64}`, time arguments of `Float64`
                                     and parameter arguments of type `Vector{Float64}` or
                                     `SciMLBase.NullParameters`.
                                     """

struct NoRecompileArgumentError <: Exception
    args
end

function Base.showerror(io::IO, e::NoRecompileArgumentError)
    println(io, NORECOMPILE_ARGUMENT_MESSAGE)
    print(io, "Attempted arguments: ")
    print(io, e.args)
end

function wrapfun_oop(ff, inputs::Tuple)
    void(ff)
end

function wrapfun_iip(ff, inputs::Tuple)
    void(ff)
end

function unwrap_fw(fw::FunctionWrapper)
    fw.obj[]
end
