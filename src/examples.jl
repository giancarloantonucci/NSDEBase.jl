"""
    Dahlquist(; u0 = 1.0, tspan = (0.0, 1.0), λ = 1.0) -> InitialValueProblem

returns an `InitialValueProblem` for the Dahlquist equation.
"""
function Dahlquist(; u0 = 1.0, tspan = (0.0, 1.0), λ = 1.0)
    f!(du, u, t) = @. du = λ * u
    return IVP(f!, u0, tspan)
end

"""
    Logistic(; u0 = 0.5, tspan = (0.0, 5.0), λ = 1.0) -> InitialValueProblem

returns an `InitialValueProblem` for the Logistic equation.
"""
function Logistic(; u0 = 0.5, tspan = (0.0, 5.0), λ = 1.0)
    f!(du, u, t) = @. du = λ * u * (1.0 - u)
    return IVP(f!, u0, tspan)
end

"""
    Riccati(; u0 = 0.0, tspan = (0.0, 5.0)) -> InitialValueProblem

returns an `InitialValueProblem` for the Riccati equation.
"""
function Riccati(; u0 = 0.0, tspan = (0.0, 5.0))
    f!(du, u, t) = @. du = u - t * u^2 + t
    return IVP(f!, u0, tspan)
end

"""
    SimplePendulum(; u0 = [0.0, π/2], tspan = (0.0, 2π), g = 1.0, L = 1.0) -> InitialValueProblem

returns an `InitialValueProblem` for the simple pendulum problem.
"""
function SimplePendulum(; u0 = [0.0, π/2], tspan = (0.0, 2π), g = 1.0, L = 1.0)
    function f!(du, u, t)
        du[1] = u[2]
        du[2] = -g/L * sin(u[1])
        return du
    end
    return IVP(f!, u0, tspan)
end

"""
    VanderPol(; u0 = [1.0; 0.0], tspan = (0.0, 5.0), μ = 1.0) -> InitialValueProblem

returns an `InitialValueProblem` for the Van der Pol equation (in first-order form).
"""
function VanderPol(; u0 = [1.0; 0.0], tspan = (0.0, 5.0), μ = 1.0)
    function f!(du, u, t)
        du[1] = u[2]
        du[2] = μ * (1.0 - u[1]^2) * u[2] - u[1]
        return du
    end
    return IVP(f!, u0, tspan)
end

"""
    Lorenz(; u0 = [2.0, 3.0, -14.0], tspan = (0.0, 10.0), σ = 10.0, β = 8/3, ρ = 28.0) -> InitialValueProblem

returns an `InitialValueProblem` for the Lorenz equations.
"""
function Lorenz(; u0 = [2.0, 3.0, -14.0], tspan = (0.0, 10.0), σ = 10.0, β = 8/3, ρ = 28.0)
    function f!(du, u, t)
        du[1] = σ * (u[2] - u[1])
        du[2] = u[1] * (ρ - u[3]) - u[2]
        du[3] = u[1] * u[2] - β * u[3]
        return du
    end
    return IVP(f!, u0, tspan)
end

"""
    Rössler(; u0 = [2.0, 0.0, 0.0], tspan = (0.0, 10.0), a = 0.2, b = 0.2, c = 5.7) -> InitialValueProblem

returns an `InitialValueProblem` for the Rössler equations.
"""
function Rössler(; u0 = [2.0, 0.0, 0.0], tspan = (0.0, 10.0), a = 0.2, b = 0.2, c = 5.7)
    function f!(du, u, t)
        du[1] = - u[2] - u[3]
        du[2] = u[1] + a * u[2]
        du[3] = b + u[3] * (u[1] - c)
        return du
    end
    return IVP(f!, u0, tspan)
end
