"""
    Dahlquist(u₀ = 1.0, tspan = (0.0, 1.0); λ = 1.0) :: InitialValueProblem
    Dahlquist(; u₀, tspan, kwargs...) :: InitialValueProblem

returns an `InitialValueProblem` for the Dahlquist equation.
"""
function Dahlquist(u₀ = 1.0, tspan = (0.0, 1.0); λ = 1.0)
    f!(du, u, t) = @. du = λ * u
    return IVP(f!, u₀, tspan)
end
Dahlquist(; u₀ = 1.0, tspan = (0.0, 1.0), kwargs...) = Dahlquist(u₀, tspan; kwargs...)

"""
    Logistic(u₀ = 0.5, tspan = (0.0, 5.0); λ = 1.0) :: InitialValueProblem
    Logistic(; u₀, tspan, kwargs...) :: InitialValueProblem

returns an `InitialValueProblem` for the Logistic equation.
"""
function Logistic(u₀ = 0.5, tspan = (0.0, 5.0); λ = 1.0)
    f!(du, u, t) = @. du = λ * u * (1.0 - u)
    return IVP(f!, u₀, tspan)
end
Logistic(; u₀ = 0.5, tspan = (0.0, 5.0), kwargs...) = Logistic(u₀, tspan; kwargs...)

"""
    Riccati(u₀ = 0.0, tspan = (0.0, 5.0)) :: InitialValueProblem
    Riccati(; u₀, tspan, kwargs...) :: InitialValueProblem

returns an `InitialValueProblem` for the Riccati equation.
"""
function Riccati(u₀ = 0.0, tspan = (0.0, 5.0))
    f!(du, u, t) = @. du = u - t * u^2 + t
    return IVP(f!, u₀, tspan)
end
Riccati(; u₀ = 0.0, tspan = (0.0, 5.0), kwargs...) = Riccati(u₀, tspan; kwargs...)

"""
    SimplePendulum(u₀ = [0.0, π/2], tspan = (0.0, 2π); g = 1.0, L = 1.0) :: InitialValueProblem
    SimplePendulum(; u₀, tspan, kwargs...) :: InitialValueProblem

returns an `InitialValueProblem` for the simple pendulum problem.
"""
function SimplePendulum(u₀ = [0.0, π/2], tspan = (0.0, 2π); g = 1.0, L = 1.0)
    function f!(du, u, t)
        du[1] = u[2]
        du[2] = -g/L * sin(u[1])
        return du
    end
    return IVP(f!, u₀, tspan)
end
SimplePendulum(; u₀ = [0.0, π/2], tspan = (0.0, 2π), kwargs...) = SimplePendulum(u₀, tspan; kwargs...)

"""
    VanderPol(u₀ = [1.0; 0.0], tspan = (0.0, 5.0); μ = 1.0) :: InitialValueProblem
    VanderPol(; u₀, tspan, kwargs...) :: InitialValueProblem

returns an `InitialValueProblem` for the Van der Pol equation (in first-order form).
"""
function VanderPol(u₀ = [1.0; 0.0], tspan = (0.0, 5.0); μ = 1.0)
    function f!(du, u, t)
        du[1] = u[2]
        du[2] = μ * (1.0 - u[1]^2) * u[2] - u[1]
        return du
    end
    return IVP(f!, u₀, tspan)
end
VanderPol(; u₀ = [1.0; 0.0], tspan = (0.0, 5.0), kwargs...) = VanderPol(u₀, tspan; kwargs...)

"""
    Lorenz(u₀ = [2.0, 3.0, -14.0], tspan = (0.0, 10.0); σ = 10.0, β = 8/3, ρ = 28.0) :: InitialValueProblem
    Lorenz(; u₀, tspan, kwargs...) :: InitialValueProblem

returns an `InitialValueProblem` for the Lorenz equations.
"""
function Lorenz(u₀ = [2.0, 3.0, -14.0], tspan = (0.0, 10.0); σ = 10.0, β = 8/3, ρ = 28.0)
    function f!(du, u, t)
        du[1] = σ * (u[2] - u[1])
        du[2] = u[1] * (ρ - u[3]) - u[2]
        du[3] = u[1] * u[2] - β * u[3]
        return du
    end
    return IVP(f!, u₀, tspan)
end
Lorenz(; u₀ = [2.0, 3.0, -14.0], tspan = (0.0, 10.0), kwargs...) = Lorenz(u₀, tspan; kwargs...)

"""
    Rössler(u₀ = [2.0, 0.0, 0.0], tspan = (0.0, 10.0); a = 0.2, b = 0.2, c = 5.7) :: InitialValueProblem
    Rössler(; u₀, tspan, kwargs...) :: InitialValueProblem

returns an `InitialValueProblem` for the Rössler equations.
"""
function Rössler(u₀ = [2.0, 0.0, 0.0], tspan = (0.0, 10.0); a = 0.2, b = 0.2, c = 5.7)
    function f!(du, u, t)
        du[1] = - u[2] - u[3]
        du[2] = u[1] + a * u[2]
        du[3] = b + u[3] * (u[1] - c)
        return du
    end
    return IVP(f!, u₀, tspan)
end
Rössler(; u₀ = [2.0, 0.0, 0.0], tspan = (0.0, 10.0), kwargs...) = Rössler(u₀, tspan; kwargs...)
