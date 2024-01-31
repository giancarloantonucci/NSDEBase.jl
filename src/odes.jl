"""
    Dahlquist(u0=0.5, tspan=(0.0, 1.0); λ=1.0)::InitialValueProblem
    Dahlquist(u0, t0, tN; kwargs...)::InitialValueProblem

returns an [`InitialValueProblem`](@ref) for the Dahlquist equation.
"""
Dahlquist(u0=0.5, tspan=(0.0, 1.0); λ=1.0) = InitialValueProblem(λ, u0, tspan)
Dahlquist(u0, t0, tN; λ) = Dahlquist(u0, (t0, tN); λ)

"""
    Logistic(u0=0.5, tspan=(0.0, 1.0); λ=1.0)::InitialValueProblem
    Logistic(u0, t0, tN; kwargs...)::InitialValueProblem

returns an [`InitialValueProblem`](@ref) for the Logistic equation.
"""
function Logistic(u0=0.5, tspan=(0.0, 1.0); λ=1.0)
    f!(du, u, t) = @. du = λ * u * (1 - u)
    return InitialValueProblem(f!, u0, tspan)
end
Logistic(u0, t0, tN; λ) = Logistic(u0, (t0, tN); λ=λ)

"""
    SimplePendulum(u0=[π/4, 0.0], tspan=(0.0, 2π/3 * √9.81))::InitialValueProblem
    SimplePendulum(u0, t0, tN; kwargs...)::InitialValueProblem

returns an [`InitialValueProblem`](@ref) for the simple pendulum problem.
"""
function SimplePendulum(u0=[π/4, 0.0], tspan=(0.0, 2π/3 * √9.81))
    function f!(du, u, t)
        du[1] = u[2]
        du[2] = - sin(u[1])
        return du
    end
    return InitialValueProblem(f!, u0, tspan)
end
SimplePendulum(u0, t0, tN) = SimplePendulum(u0, (t0, tN))

"""
    DoublePendulum(u0=[π/4, π/4, 0.0, 0.0], tspan=(0.0, 1.0); μ=1.0, λ=1.0)::InitialValueProblem
    DoublePendulum(u0, t0, tN; kwargs...)::InitialValueProblem

returns an [`InitialValueProblem`](@ref) for the double pendulum problem.
"""
function DoublePendulum(u0=[π/4, π/4, 0.0, 0.0], tspan=(0.0, 1.0); μ=1.0, λ=1.0)
    function f!(du, u, t)
        θ₁, θ₂, ω₁, ω₂ = u
        #
        M₁₁ = 1 + μ
        M₁₂ = M₂₁ = μ * λ * cos(θ₁ - θ₂)
        M₂₂ = μ * λ^2
        M = [M₁₁ M₁₂; M₂₁ M₂₂]
        #
        v₁ = - (1+μ) * sin(θ₁) - μ * λ * ω₂^2 * sin(θ₁ - θ₂)
        v₂ = - μ * λ * sin(θ₂) + μ * λ * ω₁^2 * sin(θ₁ - θ₂)
        v = [v₁; v₂]
        #
        du[1:2] .= u[3:4]
        du[3:4] .= M \ v
        return du
    end
    return InitialValueProblem(f!, u0, tspan)
end
DoublePendulum(u0, t0, tN; μ, λ) = DoublePendulum(u0, (t0, tN); μ=μ, λ=λ)

"""
    VanDerPol(u0=[1.0, 0.0], tspan=(0.0, 1.0); μ=1.0)::InitialValueProblem
    VanDerPol(u0, t0, tN; kwargs...)::InitialValueProblem

returns an [`InitialValueProblem`](@ref) for the Van der Pol equation (in first-order form).
"""
function VanDerPol(u0=[1.0, 0.0], tspan=(0.0, 1.0); μ=1.0)
    function f!(du, u, t)
        du[1] = u[2]
        du[2] = μ * (1.0 - u[1]^2) * u[2] - u[1]
        return du
    end
    return InitialValueProblem(f!, u0, tspan)
end
VanDerPol(u0, t0, tN; μ) = VanDerPol(u0, (t0, tN); μ=μ)

"""
    Rössler(u0=[2.0, 0.0, 0.0], tspan=(0.0, 1.0); α=0.2, β=0.2, γ=5.7)::InitialValueProblem
    Rössler(u0, t0, tN; kwargs...)::InitialValueProblem

returns an [`InitialValueProblem`](@ref) for the Rössler equations.
"""
function Rössler(u0=[2.0, 0.0, 0.0], tspan=(0.0, 1.0); α=0.2, β=0.2, γ=5.7)
    function f!(du, u, t)
        du[1] = - u[2] - u[3]
        du[2] = u[1] + α * u[2]
        du[3] = β + u[3] * (u[1] - γ)
        return du
    end
    return InitialValueProblem(f!, u0, tspan)
end
Rössler(u0, t0, tN; α, β, γ) = Rössler(u0, (t0, tN); α=α, β=β, γ=γ)

"""
    Lorenz(u0=[2.0, 3.0, -14.0], tspan=(0.0, 1.0); σ=10.0, β=8/3, ρ=28.0)::InitialValueProblem
    Lorenz(u0, t0, tN; kwargs...)::InitialValueProblem

returns an [`InitialValueProblem`](@ref) for the Lorenz equations.
"""
function Lorenz(u0=[2.0, 3.0, -14.0], tspan=(0.0, 1.0); σ=10.0, β=8/3, ρ=28.0)
    function f!(du, u, t)
        du[1] = σ * (u[2] - u[1])
        du[2] = u[1] * (ρ - u[3]) - u[2]
        du[3] = u[1] * u[2] - β * u[3]
        return du
    end
    return InitialValueProblem(f!, u0, tspan)
end
Lorenz(u0, t0, tN; σ, β, ρ) = Lorenz(u0, (t0, tN); σ=σ, β=β, ρ=ρ)

"""
    Lorenz96(u0=[ones(39); 1.01], tspan=(0.0, 1.0); F=8.0)::InitialValueProblem
    Lorenz96(u0, t0, tN; kwargs...)::InitialValueProblem

returns an [`InitialValueProblem`](@ref) for the Lorenz-96 equations.
"""
function Lorenz96(u0=[1.01; ones(39)], tspan=(0.0, 1.0); F=8.0)
    N = length(u0)
    if N < 4
        error("Lorenz96 requires N ≥ 4.")
    end
    function f!(du, u, t)
        for i in eachindex(du)
            # Periodic BCs at i = 1, 2, N
            if i == 1
                du[1] = (u[2] - u[N-1]) * u[N] - u[1] + F
            elseif i == 2
                du[2] = (u[3] - u[N]) * u[1] - u[2] + F
            elseif i == N
                du[N] = (u[1] - u[N-2]) * u[N-1] - u[N] + F
            else
                du[i] = (u[i+1] - u[i-2]) * u[i-1] - u[i] + F
            end
        end
        return du
    end
    return InitialValueProblem(f!, u0, tspan)
end
Lorenz96(u0, t0, tN; F) = Lorenz96(u0, (t0, tN); F=F)
