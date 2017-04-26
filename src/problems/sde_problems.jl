type SDEProblem{uType,tType,isinplace,NP,F,G,C,MM,ND} <: AbstractSDEProblem{uType,tType,isinplace,ND}
  f::F
  g::G
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

  SDEProblem{typeof(u0),promote_type(map(typeof,tspan)...),iip,typeof(noise),typeof(f),typeof(g),typeof(callback),typeof(mass_matrix),typeof(noise_rate_prototype)}(
  f,g,u0,tspan,noise,callback,mass_matrix,noise_rate_prototype)
end

# Mu' = f[1] + f[2] + ...
type SplitSDEProblem{uType,tType,isinplace,NP,F,G,C,MM,ND} <: AbstractSDEProblem{uType,tType,isinplace,ND}
  f::F
  g::G
  u0::uType
  tspan::Tuple{tType,tType}
  noise::NP
  callback::C
  mass_matrix::MM
  noise_rate_prototype::ND
end

function SplitSDEProblem(f,g,u0,tspan; iip = isinplace(f[2],3),
                         noise = nothing,
                         noise_rate_prototype = nothing,
                         callback=nothing,mass_matrix=I)
  @assert typeof(f) <: Tuple
  SplitSDEProblem{typeof(u0),promote_type(map(typeof,tspan)...),iip,typeof(noise),typeof(f),
                  typeof(g),typeof(callback),typeof(mass_matrix),
                  typeof(noise_rate_prototype)}(f,g,u0,tspan,noise,callback,mass_matrix,noise_rate_prototype)
end

# M[i]*u[i]' = f[i]
type PartitionedSDEProblem{uType,tType,isinplace,NP,F,G,C,MM,ND} <: AbstractSDEProblem{uType,tType,isinplace,ND}
  f::F
  g::G
  u0::uType
  tspan::Tuple{tType,tType}
  noise::NP
  callback::C
  mass_matrix::MM
  noise_rate_prototype::ND
end

function PartitionedSDEProblem(f,g,u0,tspan; iip = isinplace(f[1],3),
                               noise = nothing,
                               noise_rate_prototype = nothing,
                               callback=nothing,mass_matrix=nothing)
  @assert typeof(f) <: Tuple
  @assert typeof(g) <: Tuple
  @assert typeof(u0) <: Tuple
  if mass_matrix == nothing
    _mm = ((I for i in 1:length(f))...)
  end
  PartitionedSDEProblem{typeof(u0),promote_type(map(typeof,tspan)...),iip,typeof(noise),typeof(f),typeof(g),
                        typeof(callback),typeof(_mm),typeof(noise_rate_prototype)}(
                        f,g,u0,tspan,noise,callback,_mm,noise_rate_prototype)
end
