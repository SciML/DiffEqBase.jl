struct OrdinaryDiffEqTag end

const dualT = ForwardDiff.Dual{ForwardDiff.Tag{OrdinaryDiffEqTag,Float64},Float64,1}
const arglists = (Tuple{Vector{Float64},Vector{Float64},Vector{Float64},Float64},
                  Tuple{Vector{Float64},Vector{Float64},SciMLBase.NullParameters,Float64},
                  Tuple{Vector{dualT},Vector{Float64},Vector{Float64},dualT},
                  Tuple{Vector{dualT},Vector{dualT},Vector{Float64},Float64},
                  Tuple{Vector{dualT},Vector{dualT},SciMLBase.NullParameters,Float64},
                  Tuple{Vector{dualT},Vector{Float64},SciMLBase.NullParameters,dualT})
const returnlists = ntuple(x -> Nothing, length(arglists))
function void(@nospecialize(f::Function)) 
  function f2(@nospecialize(du::Vector{Float64}), @nospecialize(u::Vector{Float64}), @nospecialize(p::Vector{Float64}), @nospecialize(t::Float64))
    f(du, u, p, t)
    nothing
  end

  function f2(@nospecialize(du::Vector{Float64}), @nospecialize(u::Vector{Float64}), @nospecialize(p::SciMLBase.NullParameters), @nospecialize(t::Float64))
    f(du, u, p, t)
    nothing
  end

  function f2(@nospecialize(du::Vector{dualT}), @nospecialize(u::Vector{dualT}), @nospecialize(p::Vector{Float64}), @nospecialize(t::Float64))
    f(du, u, p, t)
    nothing
  end

  function f2(@nospecialize(du::Vector{dualT}), @nospecialize(u::Vector{dualT}), @nospecialize(p::SciMLBase.NullParameters), @nospecialize(t::Float64))
    f(du, u, p, t)
    nothing
  end
  
  function f2(@nospecialize(du::Vector{dualT}), @nospecialize(u::Vector{Float64}), @nospecialize(p::Vector{Float64}), @nospecialize(t::dualT))
    f(du, u, p, t)
    nothing
  end

  function f2(@nospecialize(du::Vector{dualT}), @nospecialize(u::Vector{Float64}), @nospecialize(p::SciMLBase.NullParameters), @nospecialize(t::dualT))
    f(du, u, p, t)
    nothing
  end
  precompile(f2, (Vector{Float64}, Vector{Float64}, Vector{Float64}, Float64))
  precompile(f2, (Vector{Float64}, Vector{Float64}, SciMLBase.NullParameters, Float64))
  precompile(f2, (Vector{dualT}, Vector{dualT}, Vector{Float64}, Float64))
  precompile(f2, (Vector{dualT}, Vector{dualT}, SciMLBase.NullParameters, Float64))
  precompile(f2, (Vector{dualT}, Vector{Float64}, Vector{Float64}, dualT))
  precompile(f2, (Vector{dualT}, Vector{Float64}, SciMLBase.NullParameters, dualT))
  f2
end

const NORECOMPILE_FUNCTION = typeof(FunctionWrappersWrappers.FunctionWrappersWrapper(void(() -> nothing), arglists, returnlists))
wrap_norecompile(f) = FunctionWrappersWrappers.FunctionWrappersWrapper(void(f), arglists, returnlists)

function ODEFunction{iip,false}(f;
  mass_matrix=I,
  analytic=nothing,
  tgrad=nothing,
  jac=nothing,
  jvp=nothing,
  vjp=nothing,
  jac_prototype=nothing,
  sparsity=jac_prototype,
  Wfact=nothing,
  Wfact_t=nothing,
  paramjac=nothing,
  syms=nothing,
  indepsym=nothing,
  observed=SciMLBase.DEFAULT_OBSERVED,
  colorvec=nothing) where {iip}

  if jac === nothing && isa(jac_prototype, AbstractDiffEqLinearOperator)
    if iip
      jac = update_coefficients! #(J,u,p,t)
    else
      jac = (u, p, t) -> update_coefficients!(deepcopy(jac_prototype), u, p, t)
    end
  end

  if jac_prototype !== nothing && colorvec === nothing && ArrayInterface.fast_matrix_colors(jac_prototype)
    _colorvec = ArrayInterface.matrix_colors(jac_prototype)
  else
    _colorvec = colorvec
  end

  ODEFunction{iip,
    NORECOMPILE_FUNCTION,typeof(mass_matrix),typeof(analytic),typeof(tgrad),typeof(jac),
    typeof(jvp),typeof(vjp),typeof(jac_prototype),typeof(sparsity),typeof(Wfact),
    typeof(Wfact_t),typeof(paramjac),typeof(syms),typeof(indepsym),
    typeof(observed),typeof(_colorvec)}(
    wrap_norecompile(f), mass_matrix, analytic, tgrad, jac,
    jvp, vjp, jac_prototype, sparsity, Wfact,
    Wfact_t, paramjac, syms, indepsym, observed, _colorvec)
end