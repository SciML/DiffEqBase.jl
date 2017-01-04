resize!(i::DEIntegrator) = error("This method has not been implemented for the integrator")
cache_iter(i::DEIntegrator) = error("This method has not been implemented for the integrator")
terminate!(i::DEIntegrator) = error("This method has not been implemented for the integrator")
get_du(i::DEIntegrator) = error("This method has not been implemented for the integrator")
get_dt(i::DEIntegrator) = error("This method has not been implemented for the integrator")
get_proposed_dt(i::DEIntegrator) = error("This method has not been implemented for the integrator")
modify_proposed_dt!(i::DEIntegrator) = error("This method has not been implemented for the integrator")
u_unmodified!(i::DEIntegrator,bool) = error("This method has not been implemented for the integrator")

@recipe function f(integrator::DEIntegrator;plot_analytic=false,vars=nothing)

  vars = interpret_vars(vars,integrator.sol)

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
      if x[j] == 0
        push!(plot_vecs[j], integrator.t)
      elseif x[j]==1 && !(typeof(integrator.u) <: AbstractArray)
        push!(plot_vecs[j], integrator.u)
      else
        push!(plot_vecs[j], integrator.u[x[j]])
      end
    end
    add_labels!(labels,x,dims,integrator.sol)
  end

  if plot_analytic
    for x in vars
      for j in 1:dims
        if x[j] == 0
          push!(plot_vecs[j], integrator.t)
        elseif x[j]==1 && !(typeof(integrator.u) <: AbstractArray)
          push!(plot_vecs[j], integrator.sol.prob.analytic(integrator.t,integrator.sol[1]))
        else
          push!(plot_vecs[j], integrator.sol.prob.analytic(integrator.t,integrator.sol[1])[x[j]])
        end
      end
      add_labels!(labels,x,dims,integrator.sol)
    end
  end

  xflip --> integrator.tdir < 0
  seriestype --> :scatter
  linewidth --> 3
  #xtickfont --> font(11)
  #ytickfont --> font(11)
  #legendfont --> font(11)
  #guidefont  --> font(11)
  label --> reshape(labels,1,length(labels))
  (plot_vecs...)
end
