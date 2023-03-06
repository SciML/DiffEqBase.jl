module DiffEqBaseMPIExt

if isdefined(Base, :get_extension)
    using DiffEqBase
    import MPI
else
    using ..DiffEqBase
    import ..MPI
end

if isdefined(MPI, :AbstractMultiRequest)
    function DiffEqBase.anyeltypedual(::Type{T},
                                      counter = 0) where {T <: MPI.AbstractMultiRequest}
        Any
    end
end

end
