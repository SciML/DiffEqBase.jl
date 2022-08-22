struct OrdinaryDiffEqTag end

const dualT = ForwardDiff.Dual{ForwardDiff.Tag{OrdinaryDiffEqTag, Float64}, Float64, 1}
const NORECOMPILE_IIP_SUPPORTED_ARGS = (Tuple{Vector{Float64}, Vector{Float64},
                                              Vector{Float64}, Float64},
                                        Tuple{Vector{Float64}, Vector{Float64},
                                              SciMLBase.NullParameters, Float64})
const iip_arglists = (Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64},
                      Tuple{Vector{Float64}, Vector{Float64}, SciMLBase.NullParameters,
                            Float64
                            },
                      Tuple{Vector{dualT}, Vector{Float64}, Vector{Float64}, dualT},
                      Tuple{Vector{dualT}, Vector{dualT}, Vector{Float64}, Float64},
                      Tuple{Vector{dualT}, Vector{dualT}, SciMLBase.NullParameters, Float64
                            },
                      Tuple{Vector{dualT}, Vector{Float64}, SciMLBase.NullParameters, dualT
                            })
const iip_returnlists = ntuple(x -> Nothing, length(iip_arglists))

struct Void{F}
    f::F
end
function (f::Void)(args...)
    f.f(args...)
    nothing
end

const oop_arglists = (Tuple{Vector{Float64}, Vector{Float64}, Float64},
                      Tuple{Vector{Float64}, SciMLBase.NullParameters, Float64},
                      Tuple{Vector{Float64}, Vector{Float64}, dualT},
                      Tuple{Vector{dualT}, Vector{Float64}, Float64},
                      Tuple{Vector{dualT}, SciMLBase.NullParameters, Float64},
                      Tuple{Vector{Float64}, SciMLBase.NullParameters, dualT})

const NORECOMPILE_OOP_SUPPORTED_ARGS = (Tuple{Vector{Float64},
                                              Vector{Float64}, Float64},
                                        Tuple{Vector{Float64},
                                              SciMLBase.NullParameters, Float64})
const oop_returnlists = (Vector{Float64}, Vector{Float64},
                         ntuple(x -> Vector{dualT}, length(oop_arglists) - 2)...)

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

function wrapfun_oop(ff, inputs::Tuple = ())
    if !isempty(inputs)
        IT = Tuple{map(typeof, inputs)...}
        if IT ∉ NORECOMPILE_OOP_SUPPORTED_ARGS
            throw(NoRecompileArgumentError(IT))
        end
    end
    FunctionWrappersWrappers.FunctionWrappersWrapper(ff, oop_arglists,
                                                     oop_returnlists)
end

function wrapfun_iip(ff, inputs::Tuple = ())
    if !isempty(inputs)
        IT = Tuple{map(typeof, inputs)...}
        if IT ∉ NORECOMPILE_IIP_SUPPORTED_ARGS
            throw(NoRecompileArgumentError(IT))
        end
    end
    FunctionWrappersWrappers.FunctionWrappersWrapper(Void(ff), iip_arglists,
                                                     iip_returnlists)
end

function unwrap_fw(fw::FunctionWrapper)
    fw.obj[]
end
