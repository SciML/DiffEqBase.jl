@require Juno="e5e0dc1b-0480-54bc-9374-aad01c23163d" begin
  # Juno Rendering
  Juno.render(i::Juno.Inline, x::DESolution) = Juno.render(i, Juno.defaultrepr(x))
  Juno.render(i::Juno.Inline, x::DEIntegrator) = Juno.render(i, Juno.defaultrepr(x))
  Juno.render(i::Juno.Inline, x::DEProblem) = Juno.render(i, Juno.defaultrepr(x))
end
