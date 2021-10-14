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
- `f :: Function` : right-hand side derivative.
- `f! :: Function` : right-hand side derivative (in-place).
- `Df :: Function` : jacobian of right-hand side derivative.
- `Df! :: Function` : jacobian of right-hand side derivative (in-place).
- `f!_or_f :: Function` : function (in-place or not) from which all other fields are automatically constructed.

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

function RightHandSideFunction(f!_or_f::Function; is_complex = false)
    if is_complex
        # Jacobian using Wirtinger derivatives
        # https://math.stackexchange.com/questions/2945446/understanding-the-chain-rule-in-the-wirtinger-calculus
        # possible alternative: use reinterpret
        if hasmethod(f!_or_f, NTuple{3, Any})
            f! = f!_or_f
            f = (u, t) -> f!_or_f(similar(u), u, t)
            Df = function (u, t)
                J_ = FiniteDifferences.jacobian(central_fdm(4, 1), u -> f(u, t), u)[1]
                n = length(u)
                P = kron(Matrix(I, n, n), 0.5 * [1.0 -1.0im; 1.0 1.0im])
                J = ((P * J_) / P)[1:2:2*n-1, 1:2:2*n-1]
                return J
            end
            Df! = function (J, du, u, t)
                J .= Df(u, t)
                return J
            end
            return RightHandSideFunction(f, f!, Df, Df!)
        # Check if f!_or_f has form f(u, t)
        elseif hasmethod(f!_or_f, NTuple{2, Any})
            f = f!_or_f
            f! = (du, u, t) -> du .= f!_or_f(u, t)
            Df = function (u, t)
                J_ = FiniteDifferences.jacobian(central_fdm(4, 1), u -> f(u, t), u)[1]
                n = length(u)
                P = kron(Matrix(I, n, n), 0.5 * [1.0 -1.0im; 1.0 1.0im])
                J = ((P * J_) / P)[1:2:2*n-1, 1:2:2*n-1]
                return J
            end
            Df! = function (J, du, u, t)
                J .= Df(u, t)
                return J
            end
            return RightHandSideFunction(f, f!, Df, Df!)
        end
    else
        # Check if f!_or_f has form f!(du, u, t)
        if hasmethod(f!_or_f, NTuple{3, Any})
            f! = f!_or_f
            f = (u, t) -> f!_or_f(similar(u), u, t)
            Df = (u, t) -> ForwardDiff.jacobian((du, u) -> f!_or_f(du, u, t), similar(u), u)
            Df! = (J, du, u, t) -> ForwardDiff.jacobian!(J, (du, u) -> f!_or_f(du, u, t), du, u)
            return RightHandSideFunction(f, f!, Df, Df!)
        # Check if f!_or_f has form f(u, t)
        elseif hasmethod(f!_or_f, NTuple{2, Any})
            f = f!_or_f
            f! = (du, u, t) -> du .= f!_or_f(u, t)
            Df = (u, t) -> ForwardDiff.jacobian(u -> f!_or_f(u, t), u)
            Df! = (J, du, u, t) -> ForwardDiff.jacobian!(J, u -> f!_or_f(u, t), u)
            return RightHandSideFunction(f, f!, Df, Df!)
        end
    end
end

@doc (@doc RightHandSideFunction) RHS(args...; kwargs...) = RightHandSideFunction(args...; kwargs...)

# ---------------------------------------------------------------------------- #
#                                   Functions                                  #
# ---------------------------------------------------------------------------- #

"""
    show(io::IO, rhs::RightHandSideFunction)

prints a full description of `rhs` and its contents to a stream `io`.
"""
Base.show(io::IO, rhs::RightHandSideFunction) = _show(io, rhs)
# function Base.show(io::IO, rhs::RightHandSideFunction)
#     print(io, "RightHandSideFunction:")
#     pad = get(io, :pad, "")
#     names = propertynames(rhs)
#     N = length(names)
#     args(n) = (n == 1 ? "(u, t)" : n == 2 ? "(du, u, t)" : n == 3 ? "(u, t)" : "(J, du, u, t)")
#     for (n, name) in enumerate(names)
#         field = getproperty(rhs, name)
#         if field !== nothing
#             print(io, "\n", pad, "   ‣ " * string(name), args(n))
#         end
#     end
# end

"""
    summary(io::IO, rhs::RightHandSideFunction)

prints a brief description of `rhs` to a stream `io`.
"""
Base.summary(io::IO, rhs::RightHandSideFunction) = _summary(io, rhs)
