### Abstract Interface

Base.length(sol::DESolution) = length(sol.u) # Must be on u for the test solutions!
Base.endof(sol::DESolution) = length(sol)
Base.getindex(sol::DESolution,i::Int) = sol.u[i]
Base.getindex(sol::DESolution,i::Int,I::Int...) = sol.u[i][I...]
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
                   plot_analytic=false,denseplot=true,plotdensity=100,vars=nothing)

  if typeof(sol) <: AbstractSDESolution
    denseplot = false
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

  function u_n(timeseries::AbstractArray, n::Int)
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

  plotx = Vector{Any}(0)
  ploty = Vector{Any}(0)
  labels = Array{String, 2}(1, length(vars)*(1+plot_analytic))
  for (i, (x, y)) in enumerate(vars)
    push!(plotx, u_n(plot_timeseries, x))
    push!(ploty, u_n(plot_timeseries, y))
    if y == 0
      ly = "t"
    else
      if has_syms(sol.prob.f)
        ly = sol.prob.f.syms[y]
      else
        ly = "u$y"
      end
    end
    if x == 0
      labels[i] = "$ly(t)"
    else
      if has_syms(sol.prob.f)
        tmp = sol.prob.f.syms[x]
        labels[i] = "($tmp, $ly)"
      else
        labels[i] = "(u$x, $ly)"
      end
    end
  end

  if plot_analytic
    for (i, (x, y)) in enumerate(vars)
      push!(plotx, u_n(plot_analytic_timeseries, x))
      push!(ploty, u_n(plot_analytic_timeseries, y))
      if y == 0
        ly = "t"
      else
        if has_syms(sol.prob.f)
          ly = string("True ",sol.prob.f.syms[y])
        else
          ly = "True u$y"
        end
      end
      if x == 0
        labels[i] = "$ly(t)"
      else
        if has_syms(sol.prob.f)
          tmp = string("True ",sol.prob.f.syms[x])
          labels[i] = "($tmp, $ly)"
        else
          labels[i] = "(True u$x, $ly)"
        end
      end
    end
  end

  seriestype --> :path
  linewidth --> 3
  #xtickfont --> font(11)
  #ytickfont --> font(11)
  #legendfont --> font(11)
  #guidefont  --> font(11)
  label --> labels
  plotx, ploty
end
