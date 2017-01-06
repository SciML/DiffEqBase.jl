### Abstract Interface

Base.length(sol::DESolution) = length(sol.u) # Must be on u for the test solutions!
Base.size(sol::DESolution,n::Int) = n==1 ? length(sol.u) : error("Only dimension 1 has a well-defined size.")
Base.endof(sol::DESolution) = length(sol)
Base.getindex(sol::DESolution,i::Int) = sol.u[i]
Base.getindex(sol::DESolution,c::AbstractArray) = [sol.u[j] for j in c]
Base.getindex(sol::DESolution,i::Int,I::Int...) = sol.u[i][I...]
Base.getindex(sol::DESolution,i::Int,::Colon,I::Int) = sol.u[i][:,I]
Base.getindex(sol::DESolution,i::Int,::Colon,I::Int...) = sol.u[i][:,I...]
Base.getindex(sol::DESolution,::Colon) = sol.u
Base.getindex(sol::DESolution,::Colon,i::Int...) = [sol.u[j][i...] for j in eachindex(sol)]
Base.getindex(sol::DESolution,c::AbstractArray,i::Int...) =  [sol.u[j][i...] for j in c]
eachindex(sol::DESolution) = eachindex(sol.t)
tuples(sol::DESolution) = tuple.(sol.t,sol.u)

function start(sol::DESolution)
  sol.tslocation = 1
  1
end

function next(sol::DESolution,state)
  state += 1
  sol.tslocation = state
  (sol,state)
end

function done(sol::DESolution,state)
  state >= length(sol)
end

function eltype(sol::DESolution)
  if typeof(sol[1]) <: AbstractArray
    return typeof(sol[1][1])
  else
    return typeof(sol[1])
  end
end

#=
function print(io::IO, sol::DESolution)
  println(io,"$(typeof(sol))")
  println(io,"u: $(sol.u)")
  println(io,"t: $(sol.t)")
  nothing
end
=#

#=
function show(io::IO,sol::DESolution)
  print(io,"$(typeof(sol))")
end
=#

@recipe function f(sol::AbstractODESolution;
                   plot_analytic=false,denseplot=true,plotdensity=10*length(sol),vars=nothing)

  if typeof(sol) <: AbstractSDESolution
    denseplot = false
  end

  vars = interpret_vars(vars,sol)

  if denseplot && sol.dense
    # Generate the points from the plot from dense function
    plott = collect(Ranges.linspace(sol.t[1],sol.t[end],plotdensity))
    plot_timeseries = sol(plott)
    if plot_analytic
      plot_analytic_timeseries = [sol.prob.analytic(t,sol.prob.u0) for t in plott]
    end
  else
    # Plot for sparse output: use the timeseries itself
    plott = sol.t
    plot_timeseries = sol.u
    if plot_analytic
      plot_analytic_timeseries = sol.u_analytic
    end
  end

  dims = length(vars[1])
  for var in vars
    @assert length(var) == dims
  end
  # Should check that all have the same dims!

  plot_vecs = []
  for i in 1:dims
    push!(plot_vecs,[])
  end
  labels = String[]# Array{String, 2}(1, length(vars)*(1+plot_analytic))
  for x in vars
    for j in 1:dims
      push!(plot_vecs[j], u_n(plot_timeseries, x[j],sol,plott,plot_timeseries))
    end
    add_labels!(labels,x,dims,sol)
  end

  if plot_analytic
    for x in vars
      for j in 1:dims
        push!(plot_vecs[j], u_n(plot_timeseries, x[j],sol,plott,plot_timeseries))
      end
      add_analytic_labels!(labels,x,dims,sol)
    end
  end

  tdir = sign(sol.t[end]-sol.t[1])
  xflip --> tdir < 0
  seriestype --> :path
  linewidth --> 3
  #xtickfont --> font(11)
  #ytickfont --> font(11)
  #legendfont --> font(11)
  #guidefont  --> font(11)
  label --> reshape(labels,1,length(labels))
  (plot_vecs...)
end

function interpret_vars(vars,sol)
  if vars != nothing && has_syms(sol.prob.f)
    # Do syms conversion
    tmp_vars = []
    for var in vars
      if typeof(var) <: Symbol
        var_int = findfirst(sol.prob.f.syms,var)
      elseif eltype(var) <: Symbol # Some kind of iterable
        var_int = map((x)->findfirst(sol.prob.f.syms,x),var)
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
    if typeof(sol[1]) <: AbstractArray
      vars = collect((0, i) for i in eachindex(sol[1]))
    else
      vars = [(0, 1)]
    end
  end
  if typeof(vars) <: Integer
    vars = [(0, vars)]
  end
  if typeof(vars) <: AbstractArray
    # If list given, its elements should be tuples, or we assume x = time
    vars = [if typeof(x) <: Tuple; x else (0, x) end for x in vars]
  end
  if typeof(vars) <: Tuple
    # If tuple given...
    if typeof(vars[1]) <: AbstractArray
      if typeof(vars[2]) <: AbstractArray
        # If both axes are lists we zip (will fail if different lengths)
        vars = collect(zip(vars[1], vars[2]))
      else
        # Just the x axis is a list
        vars = [(x, vars[2]) for x in vars[1]]
      end
    else
      if typeof(vars[2]) <: AbstractArray
        # Just the y axis is a list
        vars = [(vars[1], y) for y in vars[2]]
      else
        # Both axes are numbers
        vars = [vars]
      end
    end
  end

  # Here `vars` should be a list of tuples (x, y).
  assert(typeof(vars) <: AbstractArray)
  assert(eltype(vars) <: Tuple)
  vars
end

function add_labels!(labels,x,dims,sol)
  lys = []
  for j in 2:dims
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
  if x[1] == 0
    push!(labels,"$(lys...)(t)")
  else
    if has_syms(sol.prob.f)
      tmp = sol.prob.f.syms[x[1]]
      push!(labels,"($tmp, $(lys...))")
    else
      push!(labels,"(u$x, $(lys...))")
    end
  end
  labels
end

function add_analytic_labels!(labels,x,dims,sol)
  lys = []
  for j in 2:dims
    if x[j] == 0
      push!(lys,"t,")
    else
      if has_syms(sol.prob.f)
        push!(lys,string("True ",sol.prob.f.syms[x[j]],","))
      else
        push!(lys,"True u$(x[2]),")
      end
    end
  end
  lys[end] = lys[end][1:end-1] # Take off the last comma
  if x[1] == 0
    push!(labels,"$(lys...)(t)")
  else
    if has_syms(sol.prob.f)
      tmp = string("True ",sol.prob.f.syms[x[1]])
      push!(labels,"($tmp, $(lys...))")
    else
      push!(labels,"(True u$x, $(lys...))")
    end
  end
end

function u_n(timeseries::AbstractArray, n::Int,sol,plott,plot_timeseries)
  # Returns the nth variable from the timeseries, t if n == 0
  if n == 0
    plott
  elseif n == 1 && !(typeof(sol[1]) <: AbstractArray)
    timeseries
  else
    tmp = Vector{eltype(sol[1])}(length(plot_timeseries))
    for j in 1:length(plot_timeseries)
      tmp[j] = plot_timeseries[j][n]
    end
    tmp
  end
end
