mutable struct TimeGradientWrapper{fType,uType,P} <: Function
  f::fType
  uprev::uType
  p::P
end
(ff::TimeGradientWrapper)(t) = (du2 = similar(ff.uprev); ff.f(du2,ff.uprev,ff.p,t); du2)
(ff::TimeGradientWrapper)(du2,t) = ff.f(du2,ff.uprev,ff.p,t)

mutable struct UJacobianWrapper{fType,tType,P} <: Function
  f::fType
  t::tType
  p::P
end

(ff::UJacobianWrapper)(du1,uprev) = ff.f(du1,uprev,ff.p,ff.t)
(ff::UJacobianWrapper)(uprev) = (du1 = similar(uprev); ff.f(du1,uprev,ff.p,ff.t); du1)

mutable struct TimeDerivativeWrapper{F,uType,P} <: Function
  f::F
  u::uType
  p::P
end
(ff::TimeDerivativeWrapper)(t) = ff.f(ff.u,ff.p,t)

mutable struct UDerivativeWrapper{F,tType,P} <: Function
  f::F
  t::tType
  p::P
end
(ff::UDerivativeWrapper)(u) = ff.f(u,ff.p,ff.t)

mutable struct ParamJacobianWrapper{fType,tType,uType} <: Function
  f::fType
  t::tType
  u::uType
end

function (ff::ParamJacobianWrapper)(du1,p)
  ff.f(du1,ff.u,p,ff.t)
end

function (ff::ParamJacobianWrapper)(p)
  du1 = similar(p, size(ff.u))
  ff.f(du1,ff.u,p,ff.t)
  return du1
end
