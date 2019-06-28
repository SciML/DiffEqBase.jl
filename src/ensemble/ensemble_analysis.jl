module EnsembleAnalysis

using DiffEqBase, StaticArrays

# Getters
get_timestep(sim,i) = (getindex(sol,i) for sol in sim)
get_timepoint(sim,t) = (sol(t) for sol in sim)
function componentwise_vectors_timestep(sim,i)
  arr = [get_timestep(sim,i)...]
  if typeof(arr[1]) <: AbstractArray
   return vecarr_to_vectors(VectorOfArray(arr))
  else
   return arr
  end
 end
 function componentwise_vectors_timepoint(sim,t)
   arr = [get_timepoint(sim,t)...]
   if typeof(arr[1]) <: AbstractArray
    return vecarr_to_vectors(VectorOfArray(arr))
   else
    return arr
   end
  end

# Timestep statistics
timestep_mean(sim,i) = componentwise_mean(get_timestep(sim,i))
timestep_mean(sim,::Colon) = timeseries_steps_mean(sim)
function timestep_median(sim,i)
  arr = componentwise_vectors_timestep(sim,i)
  if typeof(first(arr)) <: AbstractArray
    return reshape([median(x) for x in arr],size(sim[1][i])...)
  else
    return median(arr)
  end
end
timestep_median(sim,::Colon) = timeseries_steps_median(sim)
function timestep_quantile(sim,q,i)
  arr = componentwise_vectors_timestep(sim,i)
  if typeof(first(arr)) <: AbstractArray
    return reshape([quantile(x,q) for x in arr],size(sim[1][i])...)
  else
    return quantile(arr,q)
  end
end
timestep_quantile(sim,q,::Colon) = timeseries_steps_quantile(sim,q)
timestep_meanvar(sim,i) = componentwise_meanvar(get_timestep(sim,i))
timestep_meanvar(sim,::Colon) = timeseries_steps_meanvar(sim)
timestep_meancov(sim,i,j) = componentwise_meancov(get_timestep(sim,i),get_timestep(sim,j))
timestep_meancov(sim,::Colon,::Colon) = timeseries_steps_meancov(sim)
timestep_meancor(sim,i,j) = componentwise_meancor(get_timestep(sim,i),get_timestep(sim,j))
timestep_meancor(sim,::Colon,::Colon) = timeseries_steps_meancor(sim)
timestep_weighted_meancov(sim,W,i,j) = componentwise_weighted_meancov(get_timestep(sim,i),get_timestep(sim,j),W)
timestep_weighted_meancov(sim,W,::Colon,::Colon) = timeseries_steps_weighted_meancov(sim,W)

function timeseries_steps_mean(sim)
  DiffEqArray([timestep_mean(sim,i) for i in 1:length(sim[1])],sim[1].t)
end
function timeseries_steps_median(sim)
  DiffEqArray([timestep_median(sim,i) for i in 1:length(sim[1])],sim[1].t)
end
function timeseries_steps_quantile(sim,q)
  DiffEqArray([timestep_quantile(sim,q,i) for i in 1:length(sim[1])],sim[1].t)
end
function timeseries_steps_meanvar(sim)
  means = typeof(sim[1][1])[]
  vars = typeof(sim[1][1])[]
  for i in 1:length(sim[1])
    m,v = timestep_meanvar(sim,i)
    push!(means,m)
    push!(vars,v)
  end
  DiffEqArray(means,sim[1].t),DiffEqArray(vars,sim[1].t)
end
function timeseries_steps_meancov(sim)
  reshape([timestep_meancov(sim,i,j) for i in 1:length(sim[1]) for j in 1:length(sim[1])],length(sim[1]),length(sim[1]))
end
function timeseries_steps_meancor(sim)
  reshape([timestep_meancor(sim,i,j) for i in 1:length(sim[1]) for j in 1:length(sim[1])],length(sim[1]),length(sim[1]))
end
function timeseries_steps_weighted_meancov(sim,W)
  reshape([timestep_meancov(sim,W,i,j) for i in 1:length(sim[1]) for j in 1:length(sim[1])],length(sim[1]),length(sim[1]))
end

timepoint_mean(sim,t) = componentwise_mean(get_timepoint(sim,t))
function timepoint_median(sim,t)
  arr = componentwise_vectors_timepoint(sim,t)
  if typeof(first(arr)) <: AbstractArray
    return reshape([median(x) for x in arr],size(sim[1][1])...)
  else
    return median(arr)
  end
end
function timepoint_quantile(sim,q,t)
  arr = componentwise_vectors_timepoint(sim,t)
  if typeof(first(arr)) <: AbstractArray
    return reshape([quantile(x,q) for x in arr],size(sim[1][1])...)
  else
    return quantile(arr,q)
  end
end
timepoint_meanvar(sim,t) = componentwise_meanvar(get_timepoint(sim,t))
timepoint_meancov(sim,t1,t2) = componentwise_meancov(get_timepoint(sim,t1),get_timepoint(sim,t2))
timepoint_meancor(sim,t1,t2) = componentwise_meancor(get_timepoint(sim,t1),get_timepoint(sim,t2))
timepoint_weighted_meancov(sim,W,t1,t2) = componentwise_weighted_meancov(get_timepoint(sim,t1),get_timepoint(sim,t2),W)

function EnsembleSummary(sim::DiffEqBase.AbstractEnsembleSolution{T,N},
                                t=sim[1].t;quantiles=[0.05,0.95]) where {T,N}
  if typeof(sim[1]) <: DiffEqBase.DESolution
    m,v = timeseries_point_meanvar(sim,t)
    qlow = timeseries_point_quantile(sim,quantiles[1],t)
    qhigh = timeseries_point_quantile(sim,quantiles[2],t)
  else
    m,v = timeseries_steps_meanvar(sim)
    qlow = timeseries_steps_quantile(sim,quantiles[1])
    qhigh = timeseries_steps_quantile(sim,quantiles[2])
  end

  num_monte = length(sim)
  EnsembleSummary{T,N,typeof(t),typeof(m),typeof(v),typeof(qlow),typeof(qhigh)}(t,m,v,qlow,qhigh,num_monte,sim.elapsedTime,sim.converged)
end

function timeseries_point_mean(sim,ts)
  DiffEqArray([timepoint_mean(sim,t) for t in ts],ts)
end
function timeseries_point_median(sim,ts)
  DiffEqArray([timepoint_median(sim,t) for t in ts],ts)
end
function timeseries_point_quantile(sim,q,ts)
  DiffEqArray([timepoint_quantile(sim,q,t) for t in ts],ts)
end
function timeseries_point_meanvar(sim,ts)
  means = typeof(sim[1][1])[]
  vars = typeof(sim[1][1])[]
  for t in ts
    m,v = timepoint_meanvar(sim,t)
    push!(means,m)
    push!(vars,v)
  end
  DiffEqArray(means,ts),DiffEqArray(vars,ts)
end
function timeseries_point_meancov(sim,ts1,ts2)
  reshape([timepoint_meancov(sim,t1,t2) for t1 in ts1 for t2 in ts2],length(ts1),length(ts2))
end
function timeseries_point_meancor(sim,ts1,ts2)
  reshape([timepoint_meancor(sim,t1,t2) for t1 in ts1 for t2 in ts2],length(ts1),length(ts2))
end
function timeseries_point_weighted_meancov(sim,W,ts1,ts2)
  reshape([timepoint_meancov(sim,W,t1,t2) for t1 in ts1 for t2 in ts2],length(ts1),length(ts2))
end

function componentwise_mean(A)
  x0 = first(A)
  n = 0
  mean = zero(x0)
  for x in A
    n += 1
    if typeof(x0) <: AbstractArray && !(typeof(x0) <: SArray)
      mean .+= x
    else
      mean += x
    end
  end
  if typeof(x0) <: AbstractArray && !(typeof(x0) <: SArray)
    mean ./= n
  else
    mean /= n
  end
  mean
end

# Welford algorithm
# https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
function componentwise_meanvar(A;bessel=true)
  x0 = first(A)
  n = 0
  mean = zero(x0)
  M2 = zero(x0)
  delta = zero(x0)
  delta2 = zero(x0)
  for x in A
    n += 1
    if typeof(x0) <: AbstractArray && !(typeof(x0) <: SArray)
      delta .= x .- mean
      mean .+= delta./n
      delta2 .= x .- mean
      M2 .+= delta.*delta2
    else
      delta = x .- mean
      mean += delta./n
      delta2 = x .- mean
      M2 += delta.*delta2
    end
  end
  if n < 2
    return NaN
  else
    if bessel
      if typeof(x0) <: AbstractArray && !(typeof(x0) <: SArray)
        M2 .= M2 ./ (n .- 1)
      else
        M2 = M2 ./ (n .- 1)
      end
    else
      if typeof(x0) <: AbstractArray && !(typeof(x0) <: SArray)
        M2 .= M2 ./ n
      else
        M2 = M2 ./ n
      end
    end
    return mean,M2
  end
end

function componentwise_meancov(A,B;bessel=true)
  x0 = first(A)
  y0 = first(B)
  n = 0
  meanx = zero(x0)
  meany = zero(y0)
  C = zero(x0)
  dx = zero(x0)
  for (x,y) in zip(A,B)
    n += 1
    if typeof(x0) <: AbstractArray && !(typeof(x0) <: SArray)
      dx .= x .- meanx
      meanx .+= dx./n
      meany .+= (y.-meany)./n
      C .+= dx .* (y .- meany)
    else
      dx = x .- meanx
      meanx += dx./n
      meany += (y.-meany)./n
      C += dx .* (y .- meany)
    end
  end
  if n < 2
    return NaN
  else
    if bessel
      if typeof(x0) <: AbstractArray && !(typeof(x0) <: SArray)
        C .= C ./ (n .- 1)
      else
        C = C ./ (n .- 1)
      end
    else
      if typeof(x0) <: AbstractArray && !(typeof(x0) <: SArray)
        C .= C ./ n
      else
        C = C ./ n
      end
    end
    return meanx,meany,C
  end
end

function componentwise_meancor(A,B;bessel=true)
  mx,my,cov = componentwise_meancov(A,B;bessel=bessel)
  mx,vx = componentwise_meanvar(A;bessel=bessel)
  my,vy = componentwise_meanvar(B;bessel=bessel)
  if typeof(vx) <: AbstractArray
    vx .= sqrt.(vx)
    vy .= sqrt.(vy)
  else
    vx = sqrt.(vx)
    vy = sqrt.(vy)
  end
  mx,my,cov./(vx.*vy)
end

function componentwise_weighted_meancov(A,B,W;weight_type=:reliability)
  x0 = first(A)
  y0 = first(B)
  w0 = first(W)
  n = 0
  meanx = zero(x0)
  meany = zero(y0)
  wsum = zero(w0)
  wsum2 = zero(w0)
  C = zero(x0)
  dx = zero(x0)
  for (x,y,w) in zip(A,B,W)
    n += 1
    if typeof(x0) <: AbstractArray && !(typeof(x0) <: SArray)
      wsum .+= w
      wsum2 .+= w.*w
      dx .= x .- meanx
      meanx .+= (w ./ wsum) .* dx
      meany .+= (w ./ wsum) .* (y .- meany)
      C .+= w .* dx .* (y .- meany)
    else
      wsum += w
      wsum2 += w.*w
      dx = x .- meanx
      meanx += (w ./ wsum) .* dx
      meany += (w ./ wsum) .* (y .- meany)
      C += w .* dx .* (y .- meany)
    end
  end
  if n < 2
    return NaN
  else
    if weight_type == :population
      if typeof(x0) <: AbstractArray && !(typeof(x0) <: SArray)
        C .= C ./ wsum
      else
        C = C ./ wsum
      end
    elseif weight_type == :reliability
      if typeof(x0) <: AbstractArray && !(typeof(x0) <: SArray)
        C .= C ./ (wsum .- wsum2 ./ wsum)
      else
        C = C ./ (wsum .- wsum2 ./ wsum)
      end
    elseif weight_type == :frequency
      if typeof(x0) <: AbstractArray && !(typeof(x0) <: SArray)
        C .= C ./ (wsum .- 1)
      else
        C = C ./ (wsum .- 1)
      end
    else
      error("The weight_type which was chosen is not allowed.")
    end
    return meanx,meany,C
  end
end

export get_timestep, get_timepoint, apply_timestep, apply_timepoint,
       componentwise_vectors_timestep, componentwise_vectors_timepoint

export componentwise_mean, componentwise_meanvar

export timestep_mean, timestep_median, timestep_quantile, timestep_meanvar,
       timestep_meancov, timestep_meancor, timestep_weighted_meancov

export timeseries_steps_mean, timeseries_steps_median, timeseries_steps_quantile,
       timeseries_steps_meanvar, timeseries_steps_meancov,
       timeseries_steps_meancor, timeseries_steps_weighted_meancov

export timepoint_mean, timepoint_median, timepoint_quantile,
       timepoint_meanvar, timepoint_meancov,
       timepoint_meancor, timepoint_weighted_meancov

export timeseries_point_mean, timeseries_point_median, timeseries_point_quantile,
       timeseries_point_meanvar, timeseries_point_meancov,
       timeseries_point_meancor, timeseries_point_weighted_meancov

end
