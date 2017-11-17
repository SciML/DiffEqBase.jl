mutable struct ParameterizedFunction{isinplace,F,P} <: ConstructedParameterizedFunction{isinplace}
  f::F
  params::P
  function ParameterizedFunction{iip}(f,p) where iip
    new{iip,typeof(f),typeof(p)}(f,p)
  end
end

function ParameterizedFunction(f,p)
  iip = numargs(f)>=4
  ParameterizedFunction{iip}(f,p)
end

(pf::ParameterizedFunction{true,F,P})(t,u,du) where {F,P} = pf.f(t,u,pf.params,du)
(pf::ParameterizedFunction{true,F,P})(t,u,params,du) where {F,P} = pf.f(t,u,params,du)
(pf::ParameterizedFunction{false,F,P})(t,u) where {F,P} = pf.f(t,u,pf.params)
(pf::ParameterizedFunction{false,F,P})(t,u,params) where {F,P} = pf.f(t,u,params)

mutable struct BVParameterizedFunction{isinplace,F,P} <: ConstructedParameterizedFunction{isinplace}
  f::F
  params::P
  function BVParameterizedFunction{iip}(f,p) where iip
    new{iip,typeof(f),typeof(p)}(f,p)
  end
end

function BVParameterizedFunction(f,p)
  iip = numargs(f)>=4
  BVParameterizedFunction{iip}(f,p)
end

(pf::BVParameterizedFunction{true,F,P})(t,u,du) where {F,P} = pf.f(t,u,pf.params,du)
(pf::BVParameterizedFunction{true,F,P})(t,u,params,du) where {F,P} = pf.f(t,u,params,du)
(pf::BVParameterizedFunction{false,F,P})(t,u) where {F,P} = pf.f(t,u,pf.params)
(pf::BVParameterizedFunction{false,F,P})(t,u,params) where {F,P} = pf.f(t,u,params)

mutable struct DAEParameterizedFunction{isinplace,F,P} <: ConstructedParameterizedFunction{isinplace}
  f::F
  params::P
  function DAEParameterizedFunction{iip}(f,p) where iip
    new{iip,typeof(f),typeof(p)}(f,p)
  end
end

function DAEParameterizedFunction(f,p)
  iip=numargs(f)>=5
  DAEParameterizedFunction{iip}(f,p)
end

(pf::DAEParameterizedFunction{true,F,P})(t,u,du,out) where {F,P} = pf.f(t,u,pf.params,du,out)
(pf::DAEParameterizedFunction{true,F,P})(t,u,params,du,out) where {F,P} = pf.f(t,u,params,du,out)
(pf::DAEParameterizedFunction{false,F,P})(t,u,du) where {F,P} = pf.f(t,u,pf.params,du)
(pf::DAEParameterizedFunction{false,F,P})(t,u,params,du) where {F,P} = pf.f(t,u,params,du)

mutable struct DDEParameterizedFunction{isinplace,F,P} <: ConstructedParameterizedFunction{isinplace}
  f::F
  params::P
  function DDEParameterizedFunction{iip}(f,p) where iip
    new{iip,typeof(f),typeof(p)}(f,p)
  end
end

function DDEParameterizedFunction(f,p)
  iip=numargs(f)>=5
  DDEParameterizedFunction{iip}(f,p)
end

(pf::DDEParameterizedFunction{true,F,P})(t,u,h,du) where {F,P} = pf.f(t,u,h,pf.params,du)
(pf::DDEParameterizedFunction{true,F,P})(t,u,h,params,du) where {F,P} = pf.f(t,u,h,params,du)
(pf::DDEParameterizedFunction{false,F,P})(t,u,h) where {F,P} = pf.f(t,u,h,pf.params)
(pf::DDEParameterizedFunction{false,F,P})(t,u,h,params) where {F,P} = pf.f(t,u,h,params)

mutable struct ProbParameterizedFunction{F,P} <: ConstructedParameterizedFunction{false}
  f::F
  params::P
end
(pf::ProbParameterizedFunction{F,P})(prob,i) where {F,P} = pf.f(prob,i,params)

mutable struct OutputParameterizedFunction{F,P} <: ConstructedParameterizedFunction{false}
  f::F
  params::P
end
(pf::OutputParameterizedFunction{F,P})(sol) where {F,P} = pf.f(sol,params)

### Interface

param_values(pf::ConstructedParameterizedFunction) = pf.params
num_params(pf::ConstructedParameterizedFunction) = length(pf.params)
function set_param_values!(pf::ConstructedParameterizedFunction,params)
  if typeof(pf.params) <: Number
    pf.params = params
  else
    pf.params .= params
  end
end
