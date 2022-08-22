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
function void(@nospecialize(f))
    function f2(@nospecialize(du::Vector{Float64}), @nospecialize(u::Vector{Float64}),
                @nospecialize(p::Vector{Float64}), @nospecialize(t::Float64))
        f(du, u, p, t)
        nothing
    end

    function f2(@nospecialize(du::Vector{Float64}), @nospecialize(u::Vector{Float64}),
                @nospecialize(p::SciMLBase.NullParameters), @nospecialize(t::Float64))
        f(du, u, p, t)
        nothing
    end

    function f2(@nospecialize(du::Vector{dualT}), @nospecialize(u::Vector{dualT}),
                @nospecialize(p::Vector{Float64}), @nospecialize(t::Float64))
        f(du, u, p, t)
        nothing
    end

    function f2(@nospecialize(du::Vector{dualT}), @nospecialize(u::Vector{dualT}),
                @nospecialize(p::SciMLBase.NullParameters), @nospecialize(t::Float64))
        f(du, u, p, t)
        nothing
    end

    function f2(@nospecialize(du::Vector{dualT}), @nospecialize(u::Vector{Float64}),
                @nospecialize(p::Vector{Float64}), @nospecialize(t::dualT))
        f(du, u, p, t)
        nothing
    end

    function f2(@nospecialize(du::Vector{dualT}), @nospecialize(u::Vector{Float64}),
                @nospecialize(p::SciMLBase.NullParameters), @nospecialize(t::dualT))
        f(du, u, p, t)
        nothing
    end
    precompile(f, (Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64))
    precompile(f, (Vector{Float64}, Vector{Float64}, SciMLBase.NullParameters, Float64))
    precompile(f, (Vector{dualT}, Vector{dualT}, Vector{Float64}, Float64))
    precompile(f, (Vector{dualT}, Vector{dualT}, SciMLBase.NullParameters, Float64))
    precompile(f, (Vector{dualT}, Vector{Float64}, Vector{Float64}, dualT))
    precompile(f, (Vector{dualT}, Vector{Float64}, SciMLBase.NullParameters, dualT))

    precompile(f2, (Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64))
    precompile(f2, (Vector{Float64}, Vector{Float64}, SciMLBase.NullParameters, Float64))
    precompile(f2, (Vector{dualT}, Vector{dualT}, Vector{Float64}, Float64))
    precompile(f2, (Vector{dualT}, Vector{dualT}, SciMLBase.NullParameters, Float64))
    precompile(f2, (Vector{dualT}, Vector{Float64}, Vector{Float64}, dualT))
    precompile(f2, (Vector{dualT}, Vector{Float64}, SciMLBase.NullParameters, dualT))
    f2
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

function typestablemapping(@nospecialize(f))
    function f2(@nospecialize(u::Vector{Float64}),
                @nospecialize(p::Vector{Float64}), @nospecialize(t::Float64))
        f(u, p, t)::Vector{Float64}
    end

    function f2(@nospecialize(u::Vector{Float64}),
                @nospecialize(p::SciMLBase.NullParameters), @nospecialize(t::Float64))
        f(u, p, t)::Vector{Float64}
    end

    function f2(@nospecialize(u::Vector{dualT}),
                @nospecialize(p::Vector{Float64}), @nospecialize(t::Float64))
        f(u, p, t)::Vector{dualT}
    end

    function f2(@nospecialize(u::Vector{dualT}),
                @nospecialize(p::SciMLBase.NullParameters), @nospecialize(t::Float64))
        f(u, p, t)::Vector{dualT}
    end

    function f2(@nospecialize(u::Vector{Float64}),
                @nospecialize(p::Vector{Float64}), @nospecialize(t::dualT))
        f(u, p, t)::Vector{dualT}
    end

    function f2(@nospecialize(u::Vector{Float64}),
                @nospecialize(p::SciMLBase.NullParameters), @nospecialize(t::dualT))
        f(u, p, t)::Vector{dualT}
    end
    precompile(f, (Vector{Float64}, Vector{Float64}, Float64))
    precompile(f, (Vector{Float64}, SciMLBase.NullParameters, Float64))
    precompile(f, (Vector{dualT}, Vector{Float64}, Float64))
    precompile(f, (Vector{dualT}, SciMLBase.NullParameters, Float64))
    precompile(f, (Vector{Float64}, Vector{Float64}, dualT))
    precompile(f, (Vector{Float64}, SciMLBase.NullParameters, dualT))

    precompile(f2, (Vector{Float64}, Vector{Float64}, Float64))
    precompile(f2, (Vector{Float64}, SciMLBase.NullParameters, Float64))
    precompile(f2, (Vector{dualT}, Vector{Float64}, Float64))
    precompile(f2, (Vector{dualT}, SciMLBase.NullParameters, Float64))
    precompile(f2, (Vector{Float64}, Vector{Float64}, dualT))
    precompile(f2, (Vector{Float64}, SciMLBase.NullParameters, dualT))
    f2
end

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
    FunctionWrappersWrappers.FunctionWrappersWrapper(typestablemapping(ff), oop_arglists,
                                                     oop_returnlists)
end

function wrapfun_iip(ff, inputs::Tuple = ())
    if !isempty(inputs)
        IT = Tuple{map(typeof, inputs)...}
        if IT ∉ NORECOMPILE_IIP_SUPPORTED_ARGS
            throw(NoRecompileArgumentError(IT))
        end
    end
    FunctionWrappersWrappers.FunctionWrappersWrapper(void(ff), iip_arglists,
                                                     iip_returnlists)
end

function unwrap_fw(fw::FunctionWrapper)
    fw.obj[]
end
