@doc raw"""
    NonlinearRightHandSide <: AbstractRightHandSide

A composite type for the right-hand side of an [`InitialValueProblem`](@ref) in the form $f(u, t)$.

# Constructors
```julia
NonlinearRightHandSide(f, f!, Df, Df!)
NonlinearRightHandSide(f!_or_f, iscomplex=false)
RightHandSide(args...; kwargs...)
RHS(args...; kwargs...)
```

# Arguments
`f :: Function` : $f$.
`f! :: Function` : $f$ (in-place).
`Df :: Function` : Jacobian of $f$.
`Df! :: Function` : Jacobian of $f$ (in-place).
"""
struct NonlinearRightHandSide{f_T<:Function, f!_T<:Function, Df_T<:Function, Df!_T<:Function} <: AbstractRightHandSide
    f::f_T
    f!::f!_T
    Df::Df_T
    Df!::Df!_T
end

function NonlinearRightHandSide(f!_or_f::Function; iscomplex::Bool=false)
    if iscomplex
        # Jacobian via Wirtinger derivatives
        if hasmethod(f!_or_f, NTuple{3, Any}) # f!_or_f like f!(du, u, t)?
            f! = f!_or_f
            f = (u, t) -> f!(similar(u), u, t)
            Df = function (u, t)
                J = FiniteDifferences.jacobian(central_fdm(4, 1), u -> f(u, t), u)[1]
                n = length(u)
                P = 0.5 * [1.0 -1.0im; 1.0 1.0im]
                Iₙ = Matrix(I, n, n)
                PP = kron(Iₙ, P)
                JJ = ((PP * J) / PP)[1:2:(2n-1), 1:2:(2n-1)]
                return JJ
            end
            Df! = function (J, du, u, t)
                J .= Df(u, t)
                return J
            end
            return NonlinearRightHandSide(f, f!, Df, Df!)
        elseif hasmethod(f!_or_f, NTuple{2, Any}) # f!_or_f like f(u, t)?
            f = f!_or_f
            f! = (du, u, t) -> du .= f(u, t)
            Df = function (u, t)
                J = FiniteDifferences.jacobian(central_fdm(4, 1), u -> f(u, t), u)[1]
                n = length(u)
                P = 0.5 * [1.0 -1.0im; 1.0 1.0im]
                Iₙ = Matrix(I, n, n)
                PP = kron(Iₙ, P)
                JJ = ((PP * J) / PP)[1:2:(2n-1), 1:2:(2n-1)]
                return JJ
            end
            Df! = function (J, du, u, t)
                J .= Df(u, t)
                return J
            end
            return NonlinearRightHandSide(f, f!, Df, Df!)
        else
            throw(ArgumentError("`NonlinearRightHandSide(f!_or_f; ...)` needs `f!_or_f` with signature `f!(du, u, t)` or `f(u, t)`."))
        end
    else
        # ToDo: NSDEFiniteDifference
        if hasmethod(f!_or_f, NTuple{3, Any}) # f!_or_f like f!(du, u, t)?
            f! = f!_or_f
            f = (u, t) -> f!(similar(u), u, t)
            Df = (u, t) -> ForwardDiff.jacobian((du, u) -> f!(du, u, t), similar(u), u)
            Df! = (J, du, u, t) -> ForwardDiff.jacobian!(J, (du, u) -> f!(du, u, t), du, u)
            return NonlinearRightHandSide(f, f!, Df, Df!)
        elseif hasmethod(f!_or_f, NTuple{2, Any}) # f!_or_f like f(u, t)?
            f = f!_or_f
            f! = (du, u, t) -> du .= f(u, t)
            Df = (u, t) -> ForwardDiff.jacobian(u -> f(u, t), u)
            Df! = (J, du, u, t) -> ForwardDiff.jacobian!(J, u -> f(u, t), u)
            return NonlinearRightHandSide(f, f!, Df, Df!)
        else
            throw(ArgumentError("`NonlinearRightHandSide(f!_or_f; ...)` needs `f!_or_f` with signature `f!(du, u, t)` or `f(u, t)`."))
        end
    end
end

@doc (@doc NonlinearRightHandSide) RightHandSide(args...; kwargs...) = NonlinearRightHandSide(args...; kwargs...)
@doc (@doc RightHandSide) RHS(args...; kwargs...) = RightHandSide(args...; kwargs...)

#----------------------------------------- METHODS -----------------------------------------

"""
    (rhs::NonlinearRightHandSide)(u, t)
    (rhs::NonlinearRightHandSide)(du, u, t)

computes the derivative `du` from the solution `u` and time `t`.
"""
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
