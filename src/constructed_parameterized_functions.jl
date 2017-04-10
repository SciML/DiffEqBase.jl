type ParameterizedFunction{isinplace,F,P} <: ConstructedParameterizedFunction{isinplace}
  f::F
  params::P
end

Base.@pure function ParameterizedFunction(f,p)
  isinplace = numargs(f)>=4
  ParameterizedFunction{isinplace,typeof(f),typeof(p)}(f,p)
end

(pf::ParameterizedFunction{true,F,P}){F,P}(t,u,du) = pf.f(t,u,pf.params,du)
(pf::ParameterizedFunction{true,F,P}){F,P}(t,u,params,du) = pf.f(t,u,params,du)
(pf::ParameterizedFunction{false,F,P}){F,P}(t,u) = pf.f(t,u,pf.params)
(pf::ParameterizedFunction{false,F,P}){F,P}(t,u,params) = pf.f(t,u,params)

type DAEParameterizedFunction{isinplace,F,P} <: ConstructedParameterizedFunction{isinplace}
  f::F
  params::P
end

Base.@pure function DAEParameterizedFunction(f,p)
  isinplace = numargs(f)>=5
  DAEParameterizedFunction{isinplace,typeof(f),typeof(p)}(f,p)
end

(pf::DAEParameterizedFunction{true,F,P}){F,P}(t,u,du,out) = pf.f(t,u,pf.params,du,out)
(pf::DAEParameterizedFunction{true,F,P}){F,P}(t,u,params,du,out) = pf.f(t,u,params,du,out)
(pf::DAEParameterizedFunction{false,F,P}){F,P}(t,u,du) = pf.f(t,u,pf.params,du)
(pf::DAEParameterizedFunction{false,F,P}){F,P}(t,u,params,du) = pf.f(t,u,params,du)

type DDEParameterizedFunction{isinplace,F,P} <: ConstructedParameterizedFunction{isinplace}
  f::F
  params::P
end

Base.@pure function DDEParameterizedFunction(f,p)
  isinplace = numargs(f)>=5
  DDEParameterizedFunction{isinplace,typeof(f),typeof(p)}(f,p)
end

(pf::DDEParameterizedFunction{true,F,P}){F,P}(t,u,h,du) = pf.f(t,u,h,pf.params,du)
(pf::DDEParameterizedFunction{true,F,P}){F,P}(t,u,h,params,du) = pf.f(t,u,h,params,du)
(pf::DDEParameterizedFunction{false,F,P}){F,P}(t,u,h) = pf.f(t,u,h,pf.params)
(pf::DDEParameterizedFunction{false,F,P}){F,P}(t,u,h,params) = pf.f(t,u,h,params)

type ProbParameterizedFunction{F,P} <: ConstructedParameterizedFunction{false}
  f::F
  params::P
end
(pf::ProbParameterizedFunction{F,P}){F,P}(prob,i) = pf.f(prob,i,params)

type OutputParameterizedFunction{F,P} <: ConstructedParameterizedFunction{false}
  f::F
  params::P
end
(pf::OutputParameterizedFunction{F,P}){F,P}(sol) = pf.f(sol,params)

### Interface

DiffEqBase.param_values(pf::ConstructedParameterizedFunction) = pf.params
DiffEqBase.num_params(pf::ConstructedParameterizedFunction) = length(pf.params)
