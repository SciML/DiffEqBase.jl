struct OrdinaryDiffEqTag end

const dualT = ForwardDiff.Dual{ForwardDiff.Tag{OrdinaryDiffEqTag, Float64}, Float64, 1}

dualgen(::Type{T}) where {T} = ForwardDiff.Dual{ForwardDiff.Tag{OrdinaryDiffEqTag, T}, T, 1}

const NORECOMPILE_IIP_SUPPORTED_ARGS = (Tuple{Vector{Float64}, Vector{Float64},
                                              Vector{Float64}, Float64},
                                        Tuple{Vector{Float64}, Vector{Float64},
                                              SciMLBase.NullParameters, Float64})

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

function wrapfun_iip(@nospecialize(ff),
                     inputs::Tuple{T1, T2, T3, T4}) where {T1, T2, T3, T4}
    T = eltype(T2)
    dualT = dualgen(T)
    dualT1 = ArrayInterfaceCore.promote_eltype(T1, dualT)
    dualT2 = ArrayInterfaceCore.promote_eltype(T2, dualT)
    dualT4 = dualgen(promote_type(T,T4))

    iip_arglists = (Tuple{T1, T2, T3, T4},
                    Tuple{dualT1, dualT2, T3, T4},
                    Tuple{dualT1, T2, T3, dualT4},
                    Tuple{dualT1, dualT2, T3, dualT4})

    iip_returnlists = ntuple(x -> Nothing, 4)

    fwt = map(iip_arglists, iip_returnlists) do A, R
        FunctionWrappersWrappers.FunctionWrappers.FunctionWrapper{R, A}(Void(ff))
    end
    FunctionWrappersWrappers.FunctionWrappersWrapper{typeof(fwt), false}(fwt)
end

const iip_arglists_default = (Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64},
                                    Float64},
                              Tuple{Vector{Float64}, Vector{Float64},
                                    SciMLBase.NullParameters,
                                    Float64
                                    },
                              Tuple{Vector{dualT}, Vector{Float64}, Vector{Float64}, dualT},
                              Tuple{Vector{dualT}, Vector{dualT}, Vector{Float64}, dualT},
                              Tuple{Vector{dualT}, Vector{dualT}, Vector{Float64}, Float64},
                              Tuple{Vector{dualT}, Vector{dualT}, SciMLBase.NullParameters,
                                    Float64
                                    },
                              Tuple{Vector{dualT}, Vector{Float64},
                                    SciMLBase.NullParameters, dualT
                                    })
const iip_returnlists_default = ntuple(x -> Nothing, length(iip_arglists_default))

function wrapfun_iip(@nospecialize(ff))
    fwt = map(iip_arglists_default, iip_returnlists_default) do A, R
        FunctionWrappersWrappers.FunctionWrappers.FunctionWrapper{R, A}(Void(ff))
    end
    FunctionWrappersWrappers.FunctionWrappersWrapper{typeof(fwt), false}(fwt)
end

function unwrap_fw(fw::FunctionWrapper)
    fw.obj[]
end
