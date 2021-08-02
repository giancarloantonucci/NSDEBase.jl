"""
    RightHandSideFunction{f_T, f!_T, Df_T, Df!_T}

returns a constructor for the right-hand side of an `InitialValueProblem`.

---

    RightHandSideFunction(f, f!, Df, Df!)
    RHS(args...; kwargs...)

returns a `RightHandSideFunction` with:
- `f   :: Function` : right-hand side derivative.
- `f!  :: Function` : right-hand side derivative (in-place).
- `Df  :: Function` : jacobian of right-hand side derivative.
- `Df! :: Function` : jacobian of right-hand side derivative (in-place).

---

    RightHandSideFunction(f!_or_f::Function)
    RHS(args...; kwargs...)

checks whether `f!_or_f` is in-place or not, and then returns a `RightHandSideFunction` with its fields constructed accordingly.
"""
struct RightHandSideFunction{f_T, f!_T, Df_T, Df!_T}
    f::f_T
    f!::f!_T
    Df::Df_T
    Df!::Df!_T
end

function RightHandSideFunction(f!_or_f::Function)
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

@doc (@doc RightHandSideFunction) RHS(args...; kwargs...) = RightHandSideFunction(args...; kwargs...)

Base.summary(io::IO, rhs::RightHandSideFunction) = print(io, "RightHandSideFunction")

function Base.show(io::IO, rhs::RightHandSideFunction)
    print(io, "RightHandSideFunction:")
    pad = get(io, :pad, "")
    names = propertynames(rhs)
    N = length(names)
    args(n) = (n == 1 ? "(u, t)" : n == 2 ? "(du, u, t)" : n == 3 ? "(u, t)" : "(J, du, u, t)")
    for (n, name) in enumerate(names)
        field = getproperty(rhs, name)
        if field !== nothing
            print(io, "\n", pad, "   ‣ " * string(name), args(n))
        end
    end
end
