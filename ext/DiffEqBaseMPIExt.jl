module DiffEqBaseMPIExt

using DiffEqBase
import MPI

if isdefined(MPI, :AbstractMultiRequest)
    function DiffEqBase.anyeltypedual(::Type{T},
            counter = 0) where {T <: MPI.AbstractMultiRequest}
        Any
    end
end

end
