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

_vec(v) = vec(v)
_vec(v::Number) = v
_vec(v::AbstractVector) = v

_reshape(v, siz) = reshape(v, siz)
_reshape(v::Number, siz) = v

macro tight_loop_macros(ex)
   :($(esc(ex)))
end

# TODO: would be good to have dtmin a function of dt
prob2dtmin(prob; use_end_time=true) = prob2dtmin(prob.tspan, oneunit(eltype(prob.tspan)), use_end_time)
function prob2dtmin(tspan, ::AbstractFloat, use_end_time)
  t1, t2 = tspan
  # handle eps(Inf) -> NaN
  t1f, t2f = map(isfinite, tspan)
  !t1f && throw(ArgumentError("t0 in the tspan `(t0, t1)` must be finite"))
  dtmin = max(eps(typeof(t1)), eps(t1))
  if use_end_time
    dtmin = t1f & t2f ? max(eps(t1), eps(t2)) : max(eps(typeof(t1)), eps(t1))
  end
  return dtmin
end
prob2dtmin(tspan, ::Integer, ::Any) = 0
# Multiplication is for putting the right units on the constant!
prob2dtmin(tspan, onet, ::Any) = onet*1//10^(10)

timedepentdtmin(integrator::DEIntegrator) = timedepentdtmin(integrator.t, integrator.opts.dtmin)
timedepentdtmin(t::AbstractFloat, dtmin) = abs(max(eps(t), dtmin))
timedepentdtmin(::Any, dtmin) = abs(dtmin)

maybe_with_logger(f, logger) = logger === nothing ? f() : Logging.with_logger(f, logger)

function default_logger(logger)
  Logging.min_enabled_level(logger) ≤ ProgressLogging.ProgressLevel && return nothing

  if Sys.iswindows() || (isdefined(Main, :IJulia) && Main.IJulia.inited)
    progresslogger = ConsoleProgressMonitor.ProgressLogger()
  else
    progresslogger = TerminalLoggers.TerminalLogger()
  end

  logger1 = LoggingExtras.EarlyFilteredLogger(progresslogger) do log
    log.level == ProgressLogging.ProgressLevel
  end
  logger2 = LoggingExtras.EarlyFilteredLogger(logger) do log
    log.level != ProgressLogging.ProgressLevel
  end

  LoggingExtras.TeeLogger(logger1, logger2)
end
