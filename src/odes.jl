"""
    Dahlquist(u0=0.5, tspan=(0.0, 1.0); λ=1.0) :: InitialValueProblem
    Dahlquist(u0, t0, tN; kwargs...) :: InitialValueProblem

returns an [`InitialValueProblem`](@ref) for the Dahlquist equation.
"""
function Dahlquist(u0=0.5, tspan=(0.0, 1.0); λ=1.0)
    f!(du, u, t) = du .= λ * u
    return IVP(f!, u0, tspan)
end
Dahlquist(u0, t0, tN; kwargs...) = Dahlquist(u0, (t0, tN); kwargs...)

"""
    Logistic(u0=0.5, tspan=(0.0, 1.0); λ=1.0) :: InitialValueProblem
    Logistic(u0, t0, tN; kwargs...) :: InitialValueProblem

returns an [`InitialValueProblem`](@ref) for the Logistic equation.
"""
function Logistic(u0=0.5, tspan=(0.0, 1.0); λ=1.0)
    f!(du, u, t) = @. du = λ * u * (1 - u)
    return IVP(f!, u0, tspan)
end
Logistic(u0, t0, tN; kwargs...) = Logistic(u0, (t0, tN); kwargs...)

"""
    SimplePendulum(u0=[π/4, 0.0], tspan=(0.0, 2π/3 * √9.81)) :: InitialValueProblem
    SimplePendulum(u0, t0, tN; kwargs...) :: InitialValueProblem

returns an [`InitialValueProblem`](@ref) for the simple pendulum problem.
"""
function SimplePendulum(u0=[π/4, 0.0], tspan=(0.0, 2π/3 * √9.81))
    function f!(du, u, t)
        du[1] = u[2]
        du[2] = - sin(u[1])
        return du
    end
    return IVP(f!, u0, tspan)
end
SimplePendulum(u0, t0, tN; kwargs...) = SimplePendulum(u0, (t0, tN); kwargs...)

"""
    DoublePendulum(u0=[π/4, π/4, 0.0, 0.0], tspan=(0.0, 1.0); μ=1.0, λ=1.0) :: InitialValueProblem
    DoublePendulum(u0, t0, tN; kwargs...) :: InitialValueProblem

returns an [`InitialValueProblem`](@ref) for the double pendulum problem.
"""
function DoublePendulum(u0=[π/4, π/4, 0.0, 0.0], tspan=(0.0, 1.0); μ=1.0, λ=1.0)
    function M(θ₁, θ₂)
      M₁₁ = 1 + μ
      M₁₂ = M₂₁ = μ * λ * cos(θ₁ - θ₂)
      M₂₂ = μ * λ^2
      return [M₁₁ M₁₂; M₂₁ M₂₂]
    end
    function v(θ₁, θ₂, ω₁, ω₂)
      v₁ = - (1+μ) * sin(θ₁) - μ * λ * ω₂^2 * sin(θ₁ - θ₂)
      v₂ = - μ * λ * sin(θ₂) + μ * λ * ω₁^2 * sin(θ₁ - θ₂)
      return [v₁; v₂]
    end
    function f!(du, u, t)
      θ₁, θ₂, ω₁, ω₂ = u
      du[1:2] = u[3:4]
      du[3:4] = M(θ₁, θ₂) \ v(θ₁, θ₂, ω₁, ω₂)
      return du
    end
    return IVP(f!, u0, tspan)
end
DoublePendulum(u0, t0, tN; kwargs...) = DoublePendulum(u0, (t0, tN); kwargs...)

"""
    VanderPol(u0=[1.0, 0.0], tspan=(0.0, 1.0); μ=1.0) :: InitialValueProblem
    VanderPol(u0, t0, tN; kwargs...) :: InitialValueProblem

returns an [`InitialValueProblem`](@ref) for the Van der Pol equation (in first-order form).
"""
function VanderPol(u0=[1.0, 0.0], tspan=(0.0, 1.0); μ=1.0)
    function f!(du, u, t)
        du[1] = u[2]
        du[2] = μ * (1.0 - u[1]^2) * u[2] - u[1]
        return du
    end
    return IVP(f!, u0, tspan)
end
VanderPol(u0, t0, tN; kwargs...) = VanderPol(u0, (t0, tN); kwargs...)

"""
    Rössler(u0=[2.0, 0.0, 0.0], tspan=(0.0, 1.0); a=0.2, b=0.2, c=5.7) :: InitialValueProblem
    Rössler(u0, t0, tN; kwargs...) :: InitialValueProblem

returns an [`InitialValueProblem`](@ref) for the Rössler equations.
"""
function Rössler(u0=[2.0, 0.0, 0.0], tspan=(0.0, 1.0); a=0.2, b=0.2, c=5.7)
    function f!(du, u, t)
        du[1] = - u[2] - u[3]
        du[2] = u[1] + a * u[2]
        du[3] = b + u[3] * (u[1] - c)
        return du
    end
    return IVP(f!, u0, tspan)
end
Rössler(u0, t0, tN; kwargs...) = Rössler(u0, (t0, tN); kwargs...)

"""
    Lorenz(u0=[2.0, 3.0, -14.0], tspan=(0.0, 1.0); σ=10.0, β=8/3, ρ=28.0) :: InitialValueProblem
    Lorenz(u0, t0, tN; kwargs...) :: InitialValueProblem

returns an [`InitialValueProblem`](@ref) for the Lorenz equations.
"""
function Lorenz(u0=[2.0, 3.0, -14.0], tspan=(0.0, 1.0); σ=10.0, β=8/3, ρ=28.0)
    function f!(du, u, t)
        du[1] = σ * (u[2] - u[1])
        du[2] = u[1] * (ρ - u[3]) - u[2]
        du[3] = u[1] * u[2] - β * u[3]
        return du
    end
    return IVP(f!, u0, tspan)
end
Lorenz(u0, t0, tN; kwargs...) = Lorenz(u0, (t0, tN); kwargs...)

"""
    Lorenz96(u0=[1.01; ones(39)], tspan=(0.0, 1.0); F=8.0) :: InitialValueProblem
    Lorenz96(u0, t0, tN; kwargs...) :: InitialValueProblem

returns an [`InitialValueProblem`](@ref) for the Lorenz-96 equations.
"""
function Lorenz96(u0=[1.01; ones(39)], tspan=(0.0, 1.0); F=8.0)
    N = length(u0)
    if N < 4
        error("Lorenz96 requires N ≥ 4.")
    end
    function f!(du, u, t)
        for i in 1:length(du)
            i == 1 ? du[1] = (u[2] - u[N - 1]) * u[N] - u[1] + F :
            i == 2 ? du[2] = (u[3] - u[N]) * u[1] - u[2] + F :
            i == N ? du[N] = (u[1] - u[N - 2]) * u[N - 1] - u[N] + F :
            du[i] = (u[i + 1] - u[i - 2]) * u[i - 1] - u[i] + F
        end
        return du
    end
    return IVP(f!, u0, tspan)
end
Lorenz96(u0, t0, tN; kwargs...) = Lorenz96(u0, (t0, tN); kwargs...)
