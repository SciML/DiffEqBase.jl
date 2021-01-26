function wrapfun_oop(ff, inputs::Tuple)
  IT = map(typeof, inputs)
  FunctionWrapper{IT[1], Tuple{IT...}}((args...)->(ff(args...)))
end

function wrapfun_iip(ff, inputs::Tuple)
  IT = map(typeof, inputs)
  FunctionWrapper{Nothing, Tuple{IT...}}((args...)->(ff(args...); nothing))
end

function unwrap_fw(fw::FunctionWrapper)
  fw.obj[]
end
