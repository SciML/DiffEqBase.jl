type RODEProblem{uType,tType,NP,F,C,MM,ND} <: AbstractRODEProblem{uType,tType,isinplace,ND}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  noise::NP
  callback::C
  mass_matrix::MM
  rand_prototype::ND
end

function RODEProblem(f,u0,tspan; iip = isinplace(f,4),
                      rand_prototype = nothing,
                      noise= nothing,
                     callback=nothing,mass_matrix=I)
  RODEProblem{typeof(u0),eltype(tspan),iip,typeof(noise),typeof(f),typeof(callback),typeof(mass_matrix),typeof(rand_prototype)}(f,u0,tspan,noise,callback,mass_matrix,rand_prototype)
end
