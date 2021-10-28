"""
    RightHandSideFunction <: AbstractRightHandSideFunction

A composite type for the right-hand side of an [`AbstractInitialValueProblem`](@ref).

# Constructors
```julia
RightHandSideFunction(f, f!, Df, Df!)
RightHandSideFunction(f!_or_f)
RHS(args...; kwargs...)
```

# Arguments
- `f::Function` : right-hand side derivative.
- `f!::Function` : right-hand side derivative (in-place).
- `Df::Function` : jacobian of right-hand side derivative.
- `Df!::Function` : jacobian of right-hand side derivative (in-place).
- `f!_or_f::Function` : function (in-place or not) from which all other fields will be constructed.

# Functions
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
struct RightHandSideFunction{f_T, f!_T, Df_T, Df!_T} <: AbstractRightHandSideFunction
    f::f_T
    f!::f!_T
    Df::Df_T
    Df!::Df!_T
end

# Check if f!_or_f is either like f!(du, u, t) or like f(u, t)
fit(f, n) = hasmethod(f, NTuple{n, Any})

function RightHandSideFunction(f!_or_f::Function; iscomplex=false)
    if iscomplex
        # Jacobian using Wirtinger derivatives:
        # https://math.stackexchange.com/questions/2945446/understanding-the-chain-rule-in-the-wirtinger-calculus
        # possible alternative: use reinterpret
        if fit(f!_or_f, 3)
            f! = f!_or_f
            f = (u, t) -> f!_or_f(similar(u), u, t)
            Df = function (u, t)
                J_ = FiniteDifferences.jacobian(central_fdm(4, 1), u -> f(u, t), u)[1]
                n = length(u)
                P̂ = 0.5 * [1.0 -1.0im; 1.0 1.0im]
                Iₙ = Matrix(I, n, n)
                P = kron(Iₙ, P̂)
                J = ((P * J_) / P)[1:2:(2 * n - 1), 1:2:(2 * n - 1)]
                return J
            end
            Df! = function (J, du, u, t)
                J .= Df(u, t)
                return J
            end
            return RightHandSideFunction(f, f!, Df, Df!)
        elseif fit(f!_or_f, 2)
            f = f!_or_f
            f! = (du, u, t) -> du .= f!_or_f(u, t)
            Df = function (u, t)
                J_ = FiniteDifferences.jacobian(central_fdm(4, 1), u -> f(u, t), u)[1]
                n = length(u)
                P̂ = 0.5 * [1.0 -1.0im; 1.0 1.0im]
                Iₙ = Matrix(I, n, n)
                P = kron(Iₙ, P̂)
                J = ((P * J_) / P)[1:2:(2 * n - 1), 1:2:(2 * n - 1)]
                return J
            end
            Df! = function (J, du, u, t)
                J .= Df(u, t)
                return J
            end
            return RightHandSideFunction(f, f!, Df, Df!)
        end
    else
        if fit(f!_or_f, 3)
            f! = f!_or_f
            f = (u, t) -> f!(similar(u), u, t)
            Df = (u, t) -> ForwardDiff.jacobian((du, u) -> f!(du, u, t), similar(u), u)
            Df! = (J, du, u, t) -> ForwardDiff.jacobian!(J, (du, u) -> f!(du, u, t), du, u)
            return RightHandSideFunction(f, f!, Df, Df!)
        elseif fit(f!_or_f, 2)
            f = f!_or_f
            f! = (du, u, t) -> du .= f(u, t)
            Df = (u, t) -> ForwardDiff.jacobian(u -> f(u, t), u)
            Df! = (J, du, u, t) -> ForwardDiff.jacobian!(J, u -> f(u, t), u)
            return RightHandSideFunction(f, f!, Df, Df!)
        end
    end
end

@doc (@doc RightHandSideFunction) function RHS(args...; kwargs...)
    return RightHandSideFunction(args...; kwargs...)
end
