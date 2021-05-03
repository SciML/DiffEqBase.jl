"""
    calculate_residuals(ũ, u₀, u₁, α, ρ, internalnorm, t)

Calculate element-wise residuals
```math
\\frac{ũ}{α+\\max{|u₀|,|u₁|}*ρ}
```
"""
@inline @muladd function calculate_residuals(ũ::Number, u₀::Number, u₁::Number,
                                    α, ρ, internalnorm,t)
    @fastmath ũ / (α + max(internalnorm(u₀,t), internalnorm(u₁,t)) * ρ)
end

@inline function calculate_residuals(ũ::Array{T}, u₀::Array{T}, u₁::Array{T}, α::T2,
                                             ρ::Real, internalnorm,t) where
                                             {T<:Number,T2<:Number}
    out = similar(ũ)
    calculate_residuals!(out, ũ, u₀, u₁, α, ρ, internalnorm,t)
    out
end

@inline function calculate_residuals(ũ, u₀, u₁, α, ρ, internalnorm,t)
  @.. calculate_residuals(ũ, u₀, u₁, α, ρ, internalnorm,t)
end

"""
    calculate_residuals(u₀, u₁, α, ρ, internalnorm, t)

Calculate element-wise residuals
```math
\\frac{u₁ - u₀}{α+\\max{|u₀|,|u₁|}*ρ}
```
"""

@inline @muladd function calculate_residuals(u₀::Number, u₁::Number,
                                    α, ρ, internalnorm,t)
    @fastmath (u₁ - u₀) / (α + max(internalnorm(u₀,t), internalnorm(u₁,t)) * ρ)
end

@inline function calculate_residuals(u₀::Array{T}, u₁::Array{T}, α::T2,
                                             ρ::Real, internalnorm,t) where
                                             {T<:Number,T2<:Number}
    out = similar(u₀)
    calculate_residuals!(out, u₀, u₁, α, ρ, internalnorm,t)
    out
end

@inline function calculate_residuals(u₀, u₁, α, ρ, internalnorm,t)
  @.. calculate_residuals(u₀, u₁, α, ρ, internalnorm,t)
end


"""
    calculate_residuals(E₁, E₂, u₀, u₁, α, ρ, δ, scalarnorm, t)

Return element-wise residuals
```math
\\frac{δ E₁ + E₂}{α+\\max{scalarnorm(u₀),scalarnorm(u₁)}*ρ}.
```
"""
@inline @muladd function calculate_residuals(E₁::Number, E₂::Number, u₀::Number, u₁::Number,
                                             α::Real, ρ::Real, δ::Number, scalarnorm, t)
    @fastmath (δ * E₁ + E₂) / (α + max(scalarnorm(u₀,t), scalarnorm(u₁,t)) * ρ)
end

@inline function calculate_residuals(E₁::Array{<:Number}, E₂::Array{<:Number},
                                     u₀::Array{<:Number}, u₁::Array{<:Number}, α::Real,
                                     ρ::Real, δ::Number, scalarnorm, t)
    out = similar(u₀)
    calculate_residuals!(out, E₁, E₂, u₀, u₁, α, ρ, δ, scalarnorm, t)
    out
end

@inline function calculate_residuals(E₁, E₂, u₀, u₁, α, ρ, δ, scalarnorm, t)
  @.. calculate_residuals(E₁, E₂, u₀, u₁, α, ρ, δ, scalarnorm, t)
end

# Inplace Versions

"""
    DiffEqBase.calculate_residuals!(out, ũ, u₀, u₁, α, ρ)

Save element-wise residuals
```math
\\frac{ũ}{α+\\max{|u₀|,|u₁|}*ρ}
```
in `out`.
"""
@inline function calculate_residuals!(out, ũ, u₀, u₁, α, ρ, internalnorm,t)
  @.. out = calculate_residuals(ũ, u₀, u₁, α, ρ, internalnorm,t)
  nothing
end


"""
    calculate_residuals!(out, u₀, u₁, α, ρ)

Save element-wise residuals
```math
\\frac{u₁ - u₀}{α+\\max{|u₀|,|u₁|}*ρ}
```
in `out`.
"""
@inline function calculate_residuals!(out, u₀, u₁, α, ρ, internalnorm,t)
  @.. out = calculate_residuals(u₀, u₁, α, ρ, internalnorm,t)
end

"""
    calculate_residuals!(out, E₁, E₂, u₀, u₁, α, ρ, δ, scalarnorm)

Calculate element-wise residuals
```math
\\frac{δ E₁ + E₂}{α+\\max{scalarnorm(u₀),scalarnorm(u₁)}*ρ}.
```
"""
@inline function calculate_residuals!(out, E₁, E₂, u₀, u₁, α, ρ, δ, scalarnorm, t)
  @.. out = calculate_residuals(E₁, E₂, u₀, u₁, α, ρ, δ, scalarnorm, t)
  out
end
