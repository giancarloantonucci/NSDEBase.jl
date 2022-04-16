@doc raw"""
    NonlinearRightHandSide <: AbstractRightHandSide

A composite type for the right-hand side of an [`InitialValueProblem`](@ref) in the form $f(u, t)$.

# Constructors
```julia
NonlinearRightHandSide(f, f!, Df, Df!)
NonlinearRightHandSide(f!_or_f, iscomplex=false)
NRHS(args...; kwargs...)
RHS(args...; kwargs...)
```
"""
struct NonlinearRightHandSide{f_T<:Function, f!_T<:Function, Df_T<:Function, Df!_T<:Function} <: AbstractRightHandSide
    f::f_T
    f!::f!_T
    Df::Df_T
    Df!::Df!_T
end
function NonlinearRightHandSide(f!_or_f::Function; iscomplex::Bool=false)
    if iscomplex
        # Jacobian via Wirtinger derivative
        # f!_or_f ?= f!(du, u, t)
        if hasmethod(f!_or_f, NTuple{3, Any})
            f! = f!_or_f
            f = (u, t) -> f!(similar(u), u, t)
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
            return NonlinearRightHandSide(f, f!, Df, Df!)
        # f!_or_f ?= f(u, t)
        elseif hasmethod(f!_or_f, NTuple{2, Any})
            f = f!_or_f
            f! = (du, u, t) -> du .= f(u, t)
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
            return NonlinearRightHandSide(f, f!, Df, Df!)
        end
    else
        # f!_or_f ?= f!(du, u, t)
        if hasmethod(f!_or_f, NTuple{3, Any})
            f! = f!_or_f
            f = (u, t) -> f!(similar(u), u, t)
            Df = (u, t) -> ForwardDiff.jacobian((du, u) -> f!(du, u, t), similar(u), u)
            Df! = (J, du, u, t) -> ForwardDiff.jacobian!(J, (du, u) -> f!(du, u, t), du, u)
            return NonlinearRightHandSide(f, f!, Df, Df!)
        # f!_or_f ?= f(u, t)
        elseif hasmethod(f!_or_f, NTuple{2, Any})
            f = f!_or_f
            f! = (du, u, t) -> du .= f(u, t)
            Df = (u, t) -> ForwardDiff.jacobian(u -> f(u, t), u)
            Df! = (J, du, u, t) -> ForwardDiff.jacobian!(J, u -> f(u, t), u)
            return NonlinearRightHandSide(f, f!, Df, Df!)
        end
    end
end
@doc (@doc NonlinearRightHandSide) NRHS(args...; kwargs...) = NonlinearRightHandSide(args...; kwargs...)
@doc (@doc NRHS) RHS(args...; kwargs...) = NRHS(args...; kwargs...)

#####
##### Methods
#####

function (rhs::NonlinearRightHandSide)(u, t)
    @↓ f = rhs
    du = f(u, t)
    return du
end
function (rhs::NonlinearRightHandSide)(du, u, t)
    @↓ f! = rhs
    # @! du = f(u, t)
    f!(du, u, t)
    return du
end
