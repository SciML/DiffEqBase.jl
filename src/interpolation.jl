immutable HermiteInterpolation{T1,T2,T3} <: AbstractDiffEqInterpolation
  t::T1
  u::T2
  du::T3
end

immutable LinearInterpolation{T1,T2} <: AbstractDiffEqInterpolation
  t::T1
  u::T2
end

interp_summary(::AbstractDiffEqInterpolation) = "Unknown"
interp_summary(::HermiteInterpolation) = "Third Order Hermite"
interp_summary(::LinearInterpolation) = "First Order Linear"

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
  if typeof(idxs) <: Number
    vals = Vector{eltype(first(timeseries))}(length(tvals))
  elseif typeof(idxs) <: AbstractVector
     vals = Vector{Vector{eltype(first(timeseries))}}(length(tvals))
  else
    vals = Vector{eltype(timeseries)}(length(tvals))
  end
  @inbounds for j in idx
    tval = tvals[j]
    i = searchsortedfirst(@view(t[i:end]),tval,rev=tdir<0)+i-1 # It's in the interval t[i-1] to t[i]
    avoid_constant_ends = deriv != Val{0} #|| typeof(tval) <: ForwardDiff.Dual
    avoid_constant_ends && i==1 && (i+=1)
    if !avoid_constant_ends && t[i] == tval
      if idxs == nothing
        vals[j] = u[i]
      else
        vals[j] = u[i][idxs]
      end
    elseif !avoid_constant_ends && t[i-1] == tval # Can happen if it's the first value!
      if idxs == nothing
        vals[j] = u[i-1]
      else
        vals[j] = u[i-1][idxs]
      end
    else
      dt = t[i] - t[i-1]
      Θ = (tval-t[i-1])/dt
      idxs_internal = idxs
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
    avoid_constant_ends = deriv != Val{0} #|| typeof(tval) <: ForwardDiff.Dual
    avoid_constant_ends && i==1 && (i+=1)
    if !avoid_constant_ends && t[i] == tval
      if idxs == nothing
        vals[j] = u[i]
      else
        vals[j] = u[i][idxs]
      end
    elseif !avoid_constant_ends && t[i-1] == tval # Can happen if it's the first value!
      if idxs == nothing
        vals[j] = u[i-1]
      else
        vals[j] = u[i-1][idxs]
      end
    else
      dt = t[i] - t[i-1]
      Θ = (tval-t[i-1])/dt
      idxs_internal = idxs
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
  avoid_constant_ends = deriv != Val{0} #|| typeof(tval) <: ForwardDiff.Dual
  avoid_constant_ends && i==1 && (i+=1)
  @inbounds if !avoid_constant_ends && t[i] == tval
    if idxs == nothing
      val = u[i]
    else
      val = u[i][idxs]
    end
  elseif !avoid_constant_ends && t[i-1] == tval # Can happen if it's the first value!
    if idxs == nothing
      val = u[i-1]
    else
      val = u[i-1][idxs]
    end
  else
    dt = t[i] - t[i-1]
    Θ = (tval-t[i-1])/dt
    idxs_internal = idxs
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
  avoid_constant_ends = deriv != Val{0} #|| typeof(tval) <: ForwardDiff.Dual
  avoid_constant_ends && i==1 && (i+=1)
  @inbounds if !avoid_constant_ends && t[i] == tval
    if idxs == nothing
      copy!(out,u[i])
    else
      copy!(out,u[i][idxs])
    end
  elseif !avoid_constant_ends && t[i-1] == tval # Can happen if it's the first value!
    if idxs == nothing
      copy!(out,u[i-1])
    else
      copy!(out,u[i-1][idxs])
    end
  else
    dt = t[i] - t[i-1]
    Θ = (tval-t[i-1])/dt
    idxs_internal = idxs
    if typeof(id) <: HermiteInterpolation
      interpolant!(out,Θ,id,dt,u[i-1],u[i],du[i-1],du[i],idxs_internal,deriv)
    else
      interpolant!(out,Θ,id,dt,u[i-1],u[i],idxs_internal,deriv)
    end
  end
end


##################### Hermite Interpolants

"""
Hairer Norsett Wanner Solving Ordinary Differential Equations I - Nonstiff Problems Page 190

Hermite Interpolation
"""
@inline function interpolant(Θ,id::HermiteInterpolation,dt,y₀,y₁,dy₀,dy₁,idxs,T::Type{Val{0}})
  if typeof(idxs) <: Void
    out = @. (1-Θ)*y₀+Θ*y₁+Θ*(Θ-1)*((1-2Θ)*(y₁-y₀)+(Θ-1)*dt*dy₀ + Θ*dt*dy₁)
  else
    out = similar(y₀,indices(idxs))
    @views @. out = (1-Θ)*y₀[idxs]+Θ*y₁[idxs]+
                    Θ*(Θ-1)*((1-2Θ)*(y₁[idxs]-y₀[idxs])+
                    (Θ-1)*dt*dy₀[idxs] + Θ*dt*dy₁[idxs])
  end
  out
end

"""
Hermite Interpolation
"""
@inline function interpolant(Θ,id::HermiteInterpolation,dt,y₀,y₁,dy₀,dy₁,idxs,T::Type{Val{1}})
  if typeof(idxs) <: Void
    out = @. dy₀ + Θ*(-4*dt*dy₀ - 2*dt*dy₁ - 6*y₀ + Θ*(3*dt*dy₀ + 3*dt*dy₁ + 6*y₀ - 6*y₁) + 6*y₁)/dt
  else
    out = similar(y₀,indices(idxs))
    @views @. out = dy₀[idxs] + Θ*(-4*dt*dy₀[idxs] -
                    2*dt*dy₁[idxs] - 6*y₀[idxs] +
                    Θ*(3*dt*dy₀[idxs] + 3*dt*dy₁[idxs] +
                    6*y₀[idxs] - 6*y₁[idxs]) + 6*y₁[idxs])/dt
  end
  out
end

"""
Hermite Interpolation
"""
@inline function interpolant(Θ,id::HermiteInterpolation,dt,y₀,y₁,dy₀,dy₁,idxs,T::Type{Val{2}})
  if typeof(idxs) <: Void
    out = @. (-4*dt*dy₀ - 2*dt*dy₁ - 6*y₀ + Θ*(6*dt*dy₀ + 6*dt*dy₁ + 12*y₀ - 12*y₁) + 6*y₁)/(dt*dt)
  else
    out = similar(y₀,indices(idxs))
    @views @. out = (-4*dt*dy₀[idxs] - 2*dt*dy₁[idxs] - 6*y₀[idxs] +
                    Θ*(6*dt*dy₀[idxs] + 6*dt*dy₁[idxs] + 12*y₀[idxs] -
                    12*y₁[idxs]) + 6*y₁[idxs])/(dt*dt)
  end
  out
end

"""
Hermite Interpolation
"""
@inline function interpolant(Θ,id::HermiteInterpolation,dt,y₀,y₁,dy₀,dy₁,idxs,T::Type{Val{3}})
  if typeof(idxs) <: Void
    out = @. (6*dt*dy₀ + 6*dt*dy₁ + 12*y₀ - 12*y₁)/(dt*dt*dt)
  else
    out = similar(y₀,indices(idxs))
    @views @. out = (6*dt*dy₀[idxs] + 6*dt*dy₁[idxs] +
                    12*y₀[idxs] - 12*y₁[idxs])/(dt*dt*dt)
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
  elseif idxs == nothing
    @. out = (1-Θ)*y₀+Θ*y₁+Θ*(Θ-1)*((1-2Θ)*(y₁-y₀)+(Θ-1)*dt*dy₀ + Θ*dt*dy₁)
  else
    @views @. out = (1-Θ)*y₀[idxs]+Θ*y₁[idxs]+Θ*(Θ-1)*((1-2Θ)*(y₁[idxs]-y₀[idxs])+(Θ-1)*dt*dy₀[idxs] + Θ*dt*dy₁[idxs])
  end
end


"""
Hermite Interpolation
"""
@inline function interpolant!(out,Θ,id::HermiteInterpolation,dt,y₀,y₁,dy₀,dy₁,idxs,T::Type{Val{1}})
  if out == nothing
    return dy₀[idxs] + Θ*(-4*dt*dy₀[idxs] - 2*dt*dy₁[idxs] - 6*y₀[idxs] + Θ*(3*dt*dy₀[idxs] + 3*dt*dy₁[idxs] + 6*y₀[idxs] - 6*y₁[idxs]) + 6*y₁[idxs])/dt
  elseif idxs == nothing
    @. out = dy₀ + Θ*(-4*dt*dy₀ - 2*dt*dy₁ - 6*y₀ + Θ*(3*dt*dy₀ + 3*dt*dy₁ + 6*y₀ - 6*y₁) + 6*y₁)/dt
  else
    @views @. out = dy₀[idxs] + Θ*(-4*dt*dy₀[idxs] - 2*dt*dy₁[idxs] - 6*y₀[idxs] + Θ*(3*dt*dy₀[idxs] + 3*dt*dy₁[idxs] + 6*y₀[idxs] - 6*y₁[idxs]) + 6*y₁[idxs])/dt
  end
end

"""
Hermite Interpolation
"""
@inline function interpolant!(out,Θ,id::HermiteInterpolation,dt,y₀,y₁,dy₀,dy₁,idxs,T::Type{Val{2}})
  if out == nothing
    return (-4*dt*dy₀[idxs] - 2*dt*dy₁[idxs] - 6*y₀[idxs] + Θ*(6*dt*dy₀[idxs] + 6*dt*dy₁[idxs] + 12*y₀[idxs] - 12*y₁[idxs]) + 6*y₁[idxs])/(dt*dt)
  elseif idxs == nothing
    @. out = (-4*dt*dy₀ - 2*dt*dy₁ - 6*y₀ + Θ*(6*dt*dy₀ + 6*dt*dy₁ + 12*y₀ - 12*y₁) + 6*y₁)/(dt*dt)
  else
    @views @. out = (-4*dt*dy₀[idxs] - 2*dt*dy₁[idxs] - 6*y₀[idxs] + Θ*(6*dt*dy₀[idxs] + 6*dt*dy₁[idxs] + 12*y₀[idxs] - 12*y₁[idxs]) + 6*y₁[idxs])/(dt*dt)
  end
end

"""
Hermite Interpolation
"""
@inline function interpolant!(out,Θ,id::HermiteInterpolation,dt,y₀,y₁,dy₀,dy₁,idxs,T::Type{Val{3}})
  if out == nothing
    return (6*dt*dy₀[idxs] + 6*dt*dy₁[idxs] + 12*y₀[idxs] - 12*y₁[idxs])/(dt*dt*dt)
  elseif idxs == nothing
    @. out = (6*dt*dy₀ + 6*dt*dy₁ + 12*y₀ - 12*y₁)/(dt*dt*dt)
  else
    @views @. out = (6*dt*dy₀[idxs] + 6*dt*dy₁[idxs] + 12*y₀[idxs] - 12*y₁[idxs])/(dt*dt*dt)
  end
end

############################### Linear Interpolants

"""
Linear Interpolation
"""
@inline function interpolant(Θ,id::LinearInterpolation,dt,y₀,y₁,idxs,T::Type{Val{0}})
  if typeof(idxs) <: Void
    out = @. (1-Θ)*y₀ + Θ*y₁
  else
    out = similar(y₀,indices(idxs))
    Θm1 = (1-Θ)
    @views @. out = Θm1*y₀[idxs] + Θ*y₁[idxs]
  end
  out
end

@inline function interpolant(Θ,id::LinearInterpolation,dt,y₀,y₁,idxs,T::Type{Val{1}})
  if typeof(idxs) <: Void
    out = @. (y₁ - y₀)/dt
  else
    out = similar(y₀,indices(idxs))
    @views @. out = (y₁[idxs] - y₀[idxs])/dt
  end
  out
end

"""
Linear Interpolation
"""
@inline function interpolant!(out,Θ,id::LinearInterpolation,dt,y₀,y₁,idxs,T::Type{Val{0}})
  Θm1 = (1-Θ)
  if out == nothing
    return Θm1*y₀[idxs] + Θ*y₁[idxs]
  elseif idxs == nothing
    @. out = Θm1*y₀ + Θ*y₁
  else
    @views @. out = Θm1*y₀[idxs] + Θ*y₁[idxs]
  end
end

"""
Linear Interpolation
"""
@inline function interpolant!(out,Θ,id::LinearInterpolation,dt,y₀,y₁,idxs,T::Type{Val{1}})
  if out == nothing
    return (y₁[idxs] - y₀[idxs])/dt
  elseif idxs == nothing
    @. out = (y₁ - y₀)/dt
  else
    @views @. out = (y₁[idxs] - y₀[idxs])/dt
  end
end
