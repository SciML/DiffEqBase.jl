module DiffEqBaseMPIExt

import DiffEqBase
isdefined(Base, :get_extension) ? (import MPI) : (import ..MPI)

if isdefined(MPI, :AbstractMultiRequest)
    function DiffEqBase.anyeltypedual(::Type{T},
        counter = 0) where {T <: MPI.AbstractMultiRequest}
        Any
    end
end

end
