immutable HermiteInterpolation{T1,T2,T3}
  t::T1
  u::T2
  du::T3
end

immutable LinearInterpolation{T1,T2}
  t::T1
  u::T2
end

(id::HermiteInterpolation)(tvals,idxs,deriv) = interpolation(tvals,id,idxs,deriv)
(id::HermiteInterpolation)(val,tvals,idxs,deriv) = interpolation!(val,tvals,id,idxs,deriv)
(id::LinearInterpolation)(tvals,idxs,deriv) = interpolation(tvals,id,idxs,deriv)
(id::LinearInterpolation)(val,tvals,idxs,deriv) = interpolation!(val,tvals,id,idxs,deriv)

@inline function interpolation(tvals,id,idxs,deriv)
  t = id.t; u = id.u
  typeof(id) <: HermiteInterpolation && (du = id.du)
  tdir = sign(t[end]-t[1])
  idx = sortperm(tvals,rev=tdir<0)
  i = 2 # Start the search thinking it's between t[1] and t[2]
  tvals[idx[end]] > t[end] && error("Solution interpolation cannot extrapolate past the final timepoint. Either solve on a longer timespan or use the local extrapolation from the integrator interface.")
  tvals[idx[1]] < t[1] && error("Solution interpolation cannot extrapolate before the first timepoint. Either start solving earlier or use the local extrapolation from the integrator interface.")
  if idxs == nothing
    if (eltype(u) <: AbstractArray) && !(eltype(u) <: Array)
      vals = Vector{Vector{eltype(first(u))}}(length(tvals))
    else
      vals = Vector{eltype(u)}(length(tvals))
    end
  elseif typeof(idxs) <: Number
    vals = Vector{eltype(first(u))}(length(tvals))
  elseif eltype(u) <: ArrayPartition
    vals = Vector{eltype(u)}(length(tvals))
  else
    vals = Vector{Vector{eltype(first(u))}}(length(tvals))
  end
  @inbounds for j in idx
    tval = tvals[j]
    i = searchsortedfirst(@view(t[i:end]),tval,rev=tdir<0)+i-1 # It's in the interval t[i-1] to t[i]
    if t[i] == tval
      if idxs == nothing
        vals[j] = u[i]
      else
        vals[j] = u[i][idxs]
      end
    elseif t[i-1] == tval # Can happen if it's the first value!
      if idxs == nothing
        vals[j] = u[i-1]
      else
        vals[j] = u[i-1][idxs]
      end
    else
      dt = t[i] - t[i-1]
      Θ = (tval-t[i-1])/dt
      if idxs == nothing && eltype(u) <: ArrayPartition
        idxs_internal = indices(u[i-1])
      elseif idxs == nothing && eltype(u) <: AbstractArray
        idxs_internal = size(u[i-1])
      else
        idxs_internal = idxs
      end
      if typeof(id) <: HermiteInterpolation
        vals[j] = interpolant(Θ,id,dt,u[i-1],u[i],du[i-1],du[i],idxs_internal,deriv)
      else
        vals[j] = interpolant(Θ,id,dt,u[i-1],u[i],idxs_internal,deriv)
      end
    end
  end
  vals
end

"""
interpolation(tvals,t,u,ks)

Get the value at tvals where the solution is known at the
times t (sorted), with values u and derivatives ks
"""
@inline function interpolation!(vals,tvals,id,idxs,deriv)
  t = id.t; u = id.u
  typeof(id) <: HermiteInterpolation && (du = id.du)
  tdir = sign(t[end]-t[1])
  idx = sortperm(tvals,rev=tdir<0)
  i = 2 # Start the search thinking it's between t[1] and t[2]
  tvals[idx[end]] > t[end] && error("Solution interpolation cannot extrapolate past the final timepoint. Either solve on a longer timespan or use the local extrapolation from the integrator interface.")
  tvals[idx[1]] < t[1] && error("Solution interpolation cannot extrapolate before the first timepoint. Either start solving earlier or use the local extrapolation from the integrator interface.")
  @inbounds for j in idx
    t = tvals[j]
    i = searchsortedfirst(@view(t[i:end]),tval,rev=tdir<0)+i-1 # It's in the interval t[i-1] to t[i]
    if t[i] == tval
      if idxs == nothing
        vals[j] = u[i]
      else
        vals[j] = u[i][idxs]
      end
    elseif t[i-1] == tval # Can happen if it's the first value!
      if idxs == nothing
        vals[j] = u[i-1]
      else
        vals[j] = u[i-1][idxs]
      end
    else
      dt = t[i] - t[i-1]
      Θ = (tval-t[i-1])/dt
      if idxs == nothing && eltype(vals) <: AbstractArray
        idxs_internal = eachindex(vals[j])
      else
        idxs_internal = idxs
      end
      if eltype(u) <: Union{AbstractArray,ArrayPartition}
        if typeof(id) <: HermiteInterpolation
          interpolant!(vals[j],Θ,id,dt,u[i-1],u[i],du[i-1],du[i],idxs_internal,deriv)
        else
          interpolant!(vals[j],Θ,id,dt,u[i-1],u[i],idxs_internal,deriv)
        end
      else
        if typeof(id) <: HermiteInterpolation
          vals[j] = interpolant(Θ,id,dt,u[i-1],u[i],du[i-1],du[i],idxs_internal,deriv)
        else
          vals[j] = interpolant(Θ,id,dt,u[i-1],u[i],idxs_internal,deriv)
        end
      end
    end
  end
end

"""
interpolation(tval::Number,t,u,ks)

Get the value at tval where the solution is known at the
times t (sorted), with values u and derivatives ks
"""
@inline function interpolation(tval::Number,id,idxs,deriv)
  t = id.t; u = id.u
  typeof(id) <: HermiteInterpolation && (du = id.du)
  tval > t[end] && error("Solution interpolation cannot extrapolate past the final timepoint. Either solve on a longer timespan or use the local extrapolation from the integrator interface.")
  tval < t[1] && error("Solution interpolation cannot extrapolate before the first timepoint. Either start solving earlier or use the local extrapolation from the integrator interface.")
  tdir = sign(t[end]-t[1])
  @inbounds i = searchsortedfirst(t,tval,rev=tdir<0) # It's in the interval t[i-1] to t[i]
  @inbounds if t[i] == tval
    if idxs == nothing
      val = u[i]
    else
      val = u[i][idxs]
    end
  elseif t[i-1] == tval # Can happen if it's the first value!
    if idxs == nothing
      val = u[i-1]
    else
      val = u[i-1][idxs]
    end
  else
    dt = t[i] - t[i-1]
    Θ = (tval-t[i-1])/dt
    if idxs == nothing && eltype(u) <: ArrayPartition
      idxs_internal = indices(u[i-1])
    elseif idxs == nothing && eltype(u) <: AbstractArray
      idxs_internal = size(u[i-1])
    else
      idxs_internal = idxs
    end
    if typeof(id) <: HermiteInterpolation
      val = interpolant(Θ,id,dt,u[i-1],u[i],du[i-1],du[i],idxs_internal,deriv)
    else
      val = interpolant(Θ,id,dt,u[i-1],u[i],idxs_internal,deriv)
    end
  end
  val
end

"""
interpolation!(out,tval::Number,t,u,ks)

Get the value at tval where the solution is known at the
times t (sorted), with values u and derivatives ks
"""
@inline function interpolation!(out,tval::Number,id,idxs,deriv)
  t = id.t; u = id.u
  typeof(id) <: HermiteInterpolation && (du = id.du)
  tval > t[end] && error("Solution interpolation cannot extrapolate past the final timepoint. Either solve on a longer timespan or use the local extrapolation from the integrator interface.")
  tval < t[1] && error("Solution interpolation cannot extrapolate before the first timepoint. Either start solving earlier or use the local extrapolation from the integrator interface.")
  tdir = sign(t[end]-t[1])
  @inbounds i = searchsortedfirst(t,tval,rev=tdir<0) # It's in the interval t[i-1] to t[i]
  @inbounds if t[i] == tval
    if idxs == nothing
      copy!(out,u[i])
    else
      copy!(out,u[i][idxs])
    end
  elseif t[i-1] == tval # Can happen if it's the first value!
    if idxs == nothing
      copy!(out,u[i-1])
    else
      copy!(out,u[i-1][idxs])
    end
  else
    dt = t[i] - t[i-1]
    Θ = (tval-t[i-1])/dt
    if idxs == nothing
      idxs_internal = eachindex(out)
    else
      idxs_internal = idxs
    end
    if typeof(id) <: HermiteInterpolation
      interpolant!(out,Θ,id,dt,u[i-1],u[i],du[i-1],du[i],idxs_internal,deriv)
    else
      interpolant!(out,Θ,id,dt,u[i-1],u[i],idxs_internal,deriv)
    end
  end
end


##################### Hermite Interpolants

#=
@inline function interpolant(Θ,id::HermiteInterpolation,dt,y₀,y₁,dy₀,dy₁,idxs,T::Type{Val{0}})
  if typeof(idxs) <: Tuple
    out = similar(y₀,idxs)
    idxs_internal=eachindex(y₀)
  else
    !(typeof(idxs) <: Number) && (out = similar(y₀,indices(idxs)))
    idxs_internal=idxs
  end
  if typeof(idxs) <: Number
    return interpolant!(nothing,Θ,dt,y₀,y₁,dy₀,dy₁,idxs_internal,T)
  else
    interpolant!(out,Θ,dt,y₀,y₁,dy₀,dy₁,idxs_internal,T)
    return out
  end
end
=#

"""
Hairer Norsett Wanner Solving Ordinary Differential Equations I - Nonstiff Problems Page 190

Hermite Interpolation
"""
@inline function interpolant(Θ,id::HermiteInterpolation,dt,y₀,y₁,dy₀,dy₁,idxs,T::Type{Val{0}}) # Default interpolant is Hermite
  if typeof(y₀) <: AbstractArray
    if typeof(idxs) <: Tuple
      out = similar(y₀,idxs)
      iter_idxs = enumerate(idxs)
    else
      out = similar(y₀,indices(idxs))
      iter_idxs = enumerate(idxs)
    end
    @inbounds for (j,i) in iter_idxs
      out[j] = (1-Θ)*y₀[i]+Θ*y₁[i]+Θ*(Θ-1)*((1-2Θ)*(y₁[i]-y₀[i])+(Θ-1)*dt*dy₀[i] + Θ*dt*dy₁[i])
    end
  else
    out = (1-Θ)*y₀+Θ*y₁+Θ*(Θ-1)*((1-2Θ)*(y₁-y₀)+(Θ-1)*dt*dy₀ + Θ*dt*dy₁)
  end
  out
end

"""
Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 190

Hermite Interpolation
"""
@inline function interpolant!(out,Θ,id::HermiteInterpolation,dt,y₀,y₁,dy₀,dy₁,idxs,T::Type{Val{0}})
  if out == nothing
    return (1-Θ)*y₀[idxs]+Θ*y₁[idxs]+Θ*(Θ-1)*((1-2Θ)*(y₁[idxs]-y₀[idxs])+(Θ-1)*dt*dy₀[idxs] + Θ*dt*dy₁[idxs])
  else
    @inbounds for (j,i) in enumerate(idxs)
      out[j] = (1-Θ)*y₀[i]+Θ*y₁[i]+Θ*(Θ-1)*((1-2Θ)*(y₁[i]-y₀[i])+(Θ-1)*dt*dy₀[i] + Θ*dt*dy₁[i])
    end
  end
end

"""
Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 190

Hermite Interpolation
"""
@inline function interpolant!(all_out::ArrayPartition,Θ,id::HermiteInterpolation,dt,all_y₀,all_y₁,all_dy₀,all_dy₁,cache,all_idxs,T::Type{Val{0}})
  for (out,y₀,y₁,idxs,dy₀,dy₁) in zip(all_out.x,all_y₀.x,all_y₁.x,all_idxs,all_dy₀.x,all_dy₁.x)
    @inbounds for (j,i) in enumerate(idxs...)
      out[j] = (1-Θ)*y₀[i]+Θ*y₁[i]+Θ*(Θ-1)*((1-2Θ)*(y₁[i]-y₀[i])+(Θ-1)*dt*dy₀ + Θ*dt*dy₁)
    end
  end
end

############################### Linear Interpolants

#=
@inline function interpolant(Θ,id::LinearInterpolation,dt,y₀,y₁,idxs,T::Type{Val{0}})
  if typeof(idxs) <: Tuple
    out = similar(y₀,idxs)
    idxs_internal=eachindex(y₀)
  else
    !(typeof(idxs) <: Number) && (out = similar(y₀,indices(idxs)))
    idxs_internal=idxs
  end
  if typeof(idxs) <: Number
    return interpolant!(nothing,Θ,dt,y₀,y₁,idxs_internal,T)
  else
    interpolant!(out,Θ,dt,y₀,y₁,idxs_internal,T)
    return out
  end
end
=#

"""
Hairer Norsett Wanner Solving Ordinary Differential Equations I - Nonstiff Problems Page 190

Linear Interpolation
"""
@inline function interpolant(Θ,id::LinearInterpolation,dt,y₀,y₁,idxs,T::Type{Val{0}}) # Default interpolant is Hermite
  if typeof(y₀) <: AbstractArray
    if typeof(idxs) <: Tuple
      out = similar(y₀,idxs)
      iter_idxs = enumerate(idxs)
    else
      out = similar(y₀,indices(idxs))
      iter_idxs = enumerate(idxs)
    end
    Θm1 = (1-Θ)
    @inbounds for (j,i) in iter_idxs
      out[j] = Θm1*y₀[i] + Θ*y₁[i]
    end
  else
    out = (1-Θ)*y₀ + Θ*y₁
  end
  out
end

"""
Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 190

Linear Interpolation
"""
@inline function interpolant!(out,Θ,id::LinearInterpolation,dt,y₀,y₁,idxs,T::Type{Val{0}})
  Θm1 = (1-Θ)
  if out == nothing
    return Θm1*y₀[idxs] + Θ*y₁[idxs]
  else
    @inbounds for (j,i) in enumerate(idxs)
      out[j] = Θm1*y₀[i] + Θ*y₁[i]
    end
  end
end

"""
Hairer Norsett Wanner Solving Ordinary Differential Euations I - Nonstiff Problems Page 190

Linear Interpolation
"""
@inline function interpolant!(all_out::ArrayPartition,Θ,id::LinearInterpolation,dt,all_y₀,all_y₁,cache,all_idxs,T::Type{Val{0}})
  Θm1 = (1-Θ)
  for (out,y₀,y₁,idxs) in zip(all_out.x,all_y₀.x,all_y₁.x,all_idxs)
    @inbounds for (j,i) in enumerate(idxs...)
      out[j] = Θm1*y₀[i] + Θ*y₁[i]
    end
  end
end
