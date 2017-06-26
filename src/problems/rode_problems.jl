type RODEProblem{uType,tType,isinplace,NP,F,C,MM,ND} <: AbstractRODEProblem{uType,tType,isinplace,ND}
  f::F
  u0::uType
  tspan::Tuple{tType,tType}
  noise::NP
  callback::C
  mass_matrix::MM
  rand_prototype::ND
  seed::UInt64
end

function RODEProblem(f,u0,tspan; iip = isinplace(f,4),
                      rand_prototype = nothing,
                      noise= nothing, seed = UInt64(0),
                     callback=nothing,mass_matrix=I)
  RODEProblem{typeof(u0),promote_type(map(typeof,tspan)...),
              iip,typeof(noise),typeof(f),typeof(callback),
              typeof(mass_matrix),typeof(rand_prototype)}(
              f,u0,tspan,noise,callback,mass_matrix,rand_prototype,seed)
end
