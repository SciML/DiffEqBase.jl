isautodifferentiable(alg::DEAlgorithm) = false
isadaptive(alg::DEAlgorithm) = true # Default to assuming adaptive, safer error("Adaptivity algorithm trait not set.")
