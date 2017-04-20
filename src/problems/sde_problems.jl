type SDEProblem{uType,tType,isinplace,NP,F,F2,C,MM,ND} <: AbstractSDEProblem{uType,tType,isinplace,ND}
  f::F
  g::F2
  u0::uType
  tspan::Tuple{tType,tType}
  noise::NP
  callback::C
  mass_matrix::MM
  noise_rate_prototype::ND
end

function SDEProblem(f,g,u0,tspan;iip = isinplace(f,3),
        noise_rate_prototype = nothing,
        noise= nothing,
        callback = nothing,mass_matrix=I)

  SDEProblem{typeof(u0),eltype(tspan),iip,typeof(noise),typeof(f),typeof(g),typeof(callback),typeof(mass_matrix),typeof(noise_rate_prototype)}(
  f,g,u0,tspan,noise,callback,mass_matrix,noise_rate_prototype)
end
