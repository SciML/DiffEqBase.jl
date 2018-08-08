### Abstract Interface

# No Time Solution : Forward to `A.u`
Base.getindex(A::AbstractNoTimeSolution,i::Int) = A.u[i]
Base.getindex(A::AbstractNoTimeSolution,I::Vararg{Int, N}) where {N} = A.u[I]
Base.setindex!(A::AbstractNoTimeSolution, v, i::Int) = (A.u[i] = v)
Base.setindex!(A::AbstractNoTimeSolution, v, I::Vararg{Int, N}) where {N} = (A.u[I] = v)
size(A::AbstractNoTimeSolution) = size(A.u)

Base.summary(A::AbstractNoTimeSolution) = string(nameof(typeof(A))," with uType ",eltype(A.u))
Base.show(io::IO, A::AbstractNoTimeSolution) = (print(io,"u: ");show(io, A.u))
Base.show(io::IO, m::MIME"text/plain", A::AbstractNoTimeSolution) = (print(io,"u: ");show(io,m,A.u))

## AbstractTimeseriesSolution Interface

Base.summary(A::AbstractTimeseriesSolution) = string(
                      TYPE_COLOR, nameof(typeof(A)),
                      NO_COLOR, " with uType ",
                      TYPE_COLOR, eltype(A.u),
                      NO_COLOR, " and tType ",
                      TYPE_COLOR, eltype(A.t), NO_COLOR)

function Base.show(io::IO, A::AbstractTimeseriesSolution)
  println(io,string("retcode: ",A.retcode))
  println(io,string("Interpolation: "),interp_summary(A.interp))
  print(io,"t: ")
  show(io, A.t)
  println(io)
  print(io,"u: ")
  show(io, A.u)
end
function Base.show(io::IO, m::MIME"text/plain", A::AbstractTimeseriesSolution)
  println(io,string("retcode: ",A.retcode))
  println(io,string("Interpolation: "),interp_summary(A.interp))
  print(io,"t: ")
  show(io,m,A.t)
  println(io)
  print(io,"u: ")
  show(io,m,A.u)
end
TreeViews.hastreeview(x::DiffEqBase.DESolution) = true
function TreeViews.treelabel(io::IO,x::DiffEqBase.DESolution,
                             mime::MIME"text/plain" = MIME"text/plain"())
  show(io,mime,Text(Base.summary(x)))
end

tuples(sol::AbstractTimeseriesSolution) = tuple.(sol.u,sol.t)

function iterate(sol::AbstractTimeseriesSolution,state=0)
  state >= length(sol) && return nothing
  state += 1
  return (solution_new_tslocation(sol,state),state)
end

DEFAULT_PLOT_FUNC(x...) = (x...,)
DEFAULT_PLOT_FUNC(x,y,z) = (x,y,z) # For v0.5.2 bug

@recipe function f(sol::AbstractTimeseriesSolution;
                   plot_analytic=false,
                   denseplot = (sol.dense ||
                                typeof(sol.prob) <: AbstractDiscreteProblem) &&
                                !(typeof(sol) <: AbstractRODESolution),
                   plotdensity = min(Int(1e5),sol.tslocation==0 ?
                                (typeof(sol.prob) <: AbstractDiscreteProblem ?
                                max(1000,100*length(sol)) :
                                max(1000,10*length(sol))) :
                                1000*sol.tslocation),
                   tspan = nothing, axis_safety = 0.1,
                   vars=nothing)

  int_vars = interpret_vars(vars,sol)

  if tspan == nothing
    if sol.tslocation == 0
      end_idx = length(sol)
    else
      end_idx = sol.tslocation
    end
    start_idx = 1
  else
    start_idx = something(findfirst(isequal(sol.t), x -> x>=tspan[end]), 0)
    end_idx = something(findfirst(isequal(sol.t), x -> x<=tspan[end]), 0)
  end

  # determine type of spacing for plott
  tscale = get(plotattributes, :xscale, :identity)
  densetspacer = if tscale in [:ln, :log10, :log2]
    (start, stop, n) -> logspace(log10(start), log10(stop), n)
  else
    linspace
  end

  if denseplot
    # Generate the points from the plot from dense function
    if tspan == nothing && !(typeof(sol) <: AbstractAnalyticalSolution)
      plott = collect(densetspacer(sol.t[start_idx],sol.t[end_idx],plotdensity))
    elseif typeof(sol) <: AbstractAnalyticalSolution
      tspan = sol.prob.tspan
      plott = collect(densetspacer(tspan[1],tspan[end],plotdensity))
    else
      plott = collect(densetspacer(tspan[1],tspan[end],plotdensity))
    end
    plot_timeseries = sol(plott)
    if plot_analytic
      if typeof(sol.prob.f) <: Tuple
        plot_analytic_timeseries = [sol.prob.f[1](Val{:analytic},t,sol.prob.u0) for t in plott]
      else
        plot_analytic_timeseries = [sol.prob.f(Val{:analytic},t,sol.prob.u0) for t in plott]
      end
    else
      plot_analytic_timeseries = nothing
    end
  else
    # Plot for sparse output: use the timeseries itself
    if sol.tslocation == 0
      plott = sol.t
      plot_timeseries = sol.u
      if plot_analytic
        plot_analytic_timeseries = sol.u_analytic
      else
        plot_analytic_timeseries = nothing
      end
    else
      if tspan == nothing
        plott = sol.t[start_idx:end_idx]
      else
        plott = collect(densetspacer(tspan[1],tspan[2],plotdensity))
      end

      plot_timeseries = sol.u[start_idx:end_idx]
      if plot_analytic
        plot_analytic_timeseries = sol.u_analytic[start_idx:end_idx]
      else
        plot_analytic_timeseries = nothing
      end
    end
  end

  dims = length(int_vars[1])
  for var in int_vars
    @assert length(var) == dims
  end
  # Should check that all have the same dims!
  plot_vecs,labels = solplot_vecs_and_labels(dims,int_vars,plot_timeseries,plott,sol,plot_analytic,plot_analytic_timeseries)

  tdir = sign(sol.t[end]-sol.t[1])
  xflip --> tdir < 0
  seriestype --> :path
  linewidth --> 3

  # Special case labels when vars = (:x,:y,:z) or (:x) or [:x,:y] ...
  if typeof(vars) <: Tuple && (typeof(vars[1]) == Symbol && typeof(vars[2]) == Symbol)
    xguide --> vars[1]
    yguide --> vars[2]
    if length(vars) > 2
      zguide --> vars[3]
    end
  end
  if getindex.(int_vars,1) == zeros(length(int_vars)) || getindex.(int_vars,2) == zeros(length(int_vars))
    xguide --> "t"
  end
  if length(int_vars[1]) >= 3 && getindex.(int_vars,3) == zeros(length(int_vars))
    yguide --> "t"
  end
  if length(int_vars[1]) >= 4 && getindex.(int_vars,4) == zeros(length(int_vars))
    zguide --> "t"
  end

  if getindex.(int_vars,2) == zeros(length(int_vars))
    if tspan == nothing
      if tdir > 0
        xlims --> (sol.t[1],sol.t[end])
      else
        xlims --> (sol.t[end],sol.t[1])
      end
    else
      xlims --> (tspan[1],tspan[end])
    end
  else
    mins = minimum(sol[int_vars[1][2],:])
    maxs = maximum(sol[int_vars[1][2],:])
    for iv in int_vars
      mins = min(mins,minimum(sol[iv[2],:]))
      maxs = max(maxs,maximum(sol[iv[2],:]))
    end
    xlims --> ((1-sign(mins)*axis_safety)*mins,(1+sign(maxs)*axis_safety)*maxs)
  end

  # Analytical solutions do not save enough information to have a good idea
  # of the axis ahead of time
  # Only set axis for animations
  if sol.tslocation != 0 && !(typeof(sol) <: AbstractAnalyticalSolution)
    if all(getindex.(int_vars,1) .== DiffEqBase.DEFAULT_PLOT_FUNC)
      mins = minimum(sol[int_vars[1][3],:])
      maxs = maximum(sol[int_vars[1][3],:])
      for iv in int_vars
        mins = min(mins,minimum(sol[iv[3],:]))
        maxs = max(maxs,maximum(sol[iv[3],:]))
      end
      ylims --> ((1-sign(mins)*axis_safety)*mins,(1+sign(maxs)*axis_safety)*maxs)

      if length(int_vars[1]) >= 4
        mins = minimum(sol[int_vars[1][4],:])
        maxs = maximum(sol[int_vars[1][4],:])
        for iv in int_vars
          mins = min(mins,minimum(sol[iv[4],:]))
          maxs = max(mins,maximum(sol[iv[4],:]))
        end
        zlims --> ((1-sign(mins)*axis_safety)*mins,(1+sign(maxs)*axis_safety)*maxs)
      end
    end
  end

  label --> reshape(labels,1,length(labels))
  (plot_vecs...,)
end

function interpret_vars(vars,sol)
  if vars != nothing && has_syms(sol.prob.f)
    # Do syms conversion
    tmp_vars = []
    for var in vars
      if typeof(var) <: Symbol
        var_int = something(findfirst(isequal(var), sol.prob.f.syms), 0)
      elseif typeof(var) <: Union{Tuple,AbstractArray} #eltype(var) <: Symbol # Some kind of iterable
        tmp = []
        for x in var
          if typeof(x) <: Symbol
            push!(tmp,something(findfirst(isequal(x), sol.prob.f.syms), 0))
          else
            push!(tmp,x)
          end
        end
        if typeof(var) <: Tuple
          var_int = tuple(tmp...)
        else
          var_int = tmp
        end
      else
        var_int = var
      end
      push!(tmp_vars,var_int)
    end
    if typeof(vars) <: Tuple
      vars = tuple(tmp_vars...)
    else
      vars = tmp_vars
    end
  end

  if vars == nothing
    # Default: plot all timeseries
    if typeof(sol[1]) <: Union{Tuple,AbstractArray}
      vars = collect((DEFAULT_PLOT_FUNC,0, i) for i in plot_indices(sol[1]))
    else
      vars = [(DEFAULT_PLOT_FUNC,0, 1)]
    end
  end

  if typeof(vars) <: Integer
    vars = [(DEFAULT_PLOT_FUNC,0, vars)]
  end

  if typeof(vars) <: AbstractArray
    # If list given, its elements should be tuples, or we assume x = time
    tmp = Tuple[]
    for x in vars
      if typeof(x) <: Tuple
        if typeof(x[1]) <: Int
          push!(tmp,tuple(DEFAULT_PLOT_FUNC,x...))
        else
          push!(tmp,x)
        end
      else
        push!(tmp,(DEFAULT_PLOT_FUNC,0, x))
      end
    end
    vars = tmp
  end

  if typeof(vars) <: Tuple
    # If tuple given...
    if typeof(vars[end-1]) <: AbstractArray
      if typeof(vars[end]) <: AbstractArray
        # If both axes are lists we zip (will fail if different lengths)
        vars = collect(zip([DEFAULT_PLOT_FUNC for i in eachindex(vars[end-1])],vars[end-1], vars[end]))
      else
        # Just the x axis is a list
        vars = [(DEFAULT_PLOT_FUNC,x, vars[end]) for x in vars[end-1]]
      end
    else
      if typeof(vars[2]) <: AbstractArray
        # Just the y axis is a list
        vars = [(DEFAULT_PLOT_FUNC,vars[end-1], y) for y in vars[end]]
      else
        # Both axes are numbers
        if typeof(vars[1]) <: Int
          vars = [tuple(DEFAULT_PLOT_FUNC,vars...)]
        else
          vars = [vars]
        end
      end
    end
  end

  # Here `vars` should be a list of tuples (x, y).
  @assert(typeof(vars) <: AbstractArray)
  @assert(eltype(vars) <: Tuple)
  vars
end

function add_labels!(labels,x,dims,sol)
  lys = []
  for j in 3:dims
    if x[j] == 0
      push!(lys,"t,")
    else
      if has_syms(sol.prob.f)
        push!(lys,"$(sol.prob.f.syms[x[j]]),")
      else
        push!(lys,"u$(x[j]),")
      end
    end
  end
  lys[end] = lys[end][1:end-1] # Take off the last comma
  if x[2] == 0 && dims == 3
    tmp_lab = "$(lys...)(t)"
  else
    if has_syms(sol.prob.f) && x[2] != 0
      tmp = sol.prob.f.syms[x[2]]
      tmp_lab = "($tmp,$(lys...))"
    else
      if x[2] == 0
        tmp_lab = "(t,$(lys...))"
      else
        tmp_lab = "(u$(x[2]),$(lys...))"
      end
    end
  end
  if x[1] != DEFAULT_PLOT_FUNC
    push!(labels,"$(x[1])$(tmp_lab)")
  else
    push!(labels,tmp_lab)
  end
  labels
end

function add_analytic_labels!(labels,x,dims,sol)
  lys = []
  for j in 3:dims
    if x[j] == 0 && dims == 3
      push!(lys,"t,")
    else
      if has_syms(sol.prob.f)
        push!(lys,string("True ",sol.prob.f.syms[x[j]],","))
      else
        push!(lys,"True u$(x[j]),")
      end
    end
  end
  lys[end] = lys[end][1:end-1] # Take off the last comma
  if x[2] == 0
    tmp_lab = "$(lys...)(t)"
  else
    if has_syms(sol.prob.f)
      tmp = string("True ",sol.prob.f.syms[x[2]])
      tmp_lab = "($tmp,$(lys...))"
    else
      tmp_lab = "(True u$(x[2]),$(lys...))"
    end
  end
  if x[1] != DEFAULT_PLOT_FUNC
    push!(labels,"$(x[1])$(tmp_lab)")
  else
    push!(labels,tmp_lab)
  end
end

function u_n(timeseries::AbstractArray, n::Int,sol,plott,plot_timeseries)
  # Returns the nth variable from the timeseries, t if n == 0
  if n == 0
    return plott
  elseif n == 1 && !(typeof(sol[1]) <: Union{AbstractArray,ArrayPartition})
    return timeseries
  else
    tmp = Vector{eltype(sol[1])}(length(plot_timeseries))
    for j in 1:length(plot_timeseries)
      tmp[j] = plot_timeseries[j][n]
    end
    return tmp
  end
end

function solplot_vecs_and_labels(dims,vars,plot_timeseries,plott,sol,plot_analytic,plot_analytic_timeseries)
  plot_vecs = []
  labels = String[]
  for x in vars
    tmp = []
    for j in 2:dims
      push!(tmp, u_n(plot_timeseries, x[j],sol,plott,plot_timeseries))
    end

    f = x[1]

    tmp = f.(tmp...)

    tmp = tuple((getindex.(tmp,i) for i in eachindex(tmp[1]))...)
    for i in eachindex(tmp)
      if length(plot_vecs) < i
        push!(plot_vecs,[])
      end
      push!(plot_vecs[i],tmp[i])
    end
    add_labels!(labels,x,dims,sol)
  end



  if plot_analytic
    analytic_plot_vecs = []
    for x in vars
      tmp = []
      for j in 2:dims
        push!(tmp, u_n(plot_analytic_timeseries, x[j],sol,plott,plot_analytic_timeseries))
      end
      f = x[1]
      tmp = f.(tmp...)
      tmp = tuple((getindex.(tmp,i) for i in eachindex(tmp[1]))...)
      for i in eachindex(tmp)
        push!(plot_vecs[i],tmp[i])
      end
      add_analytic_labels!(labels,x,dims,sol)
    end
  end
  plot_vecs = [hcat(x...) for x in plot_vecs]
  plot_vecs,labels
end

plot_indices(A::AbstractArray) = eachindex(A)
plot_indices(A::ArrayPartition) = eachindex(A)
