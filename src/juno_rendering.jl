@require Juno begin
  # Juno Rendering
  Juno.render(i::Juno.Inline, x::DESolution) = Juno.render(i, Juno.defaultrepr(x))
  Juno.render(i::Juno.Inline, x::DEIntegrator) = Juno.render(i, Juno.defaultrepr(x))
  Juno.render(i::Juno.Inline, x::DEProblem) = Juno.render(i, Juno.defaultrepr(x))
end
