@doc """
    RightHandSideFunction(f, f!, Df, Df!) -> RightHandSideFunction
    RightHandSideFunction(_f) -> RightHandSideFunction
    RHS(args...; kwargs...) -> RightHandSideFunction

returns a constructor for the right-hand side derivative of an `InitialValueProblem` with
- `f` : right-hand side derivative.
- `f!` : like `f` but in-place.
- `Df` : jacobian of `f`.
- `Df!` : like `Df` but in-place.
"""
struct RightHandSideFunction{f_T, f!_T, Df_T, Df!_T}
    f::f_T
    f!::f!_T
    Df::Df_T
    Df!::Df!_T
end

function RightHandSideFunction(_f::Function)
    if hasmethod(_f, NTuple{2, Any})
        f = _f
        f! = (du, u, t) -> du .= _f(u, t)
        Df = (u, t) -> jacobian(u -> _f(u, t), u)
        Df! = (J, du, u, t) -> jacobian!(J, u -> _f(u, t), u)
        return RightHandSideFunction(f, f!, Df, Df!)
    elseif hasmethod(_f, NTuple{3, Any})
        f! = _f
        f = (u, t) -> _f(similar(u), u, t)
        Df = (u, t) -> jacobian((du, u) -> _f(du, u, t), similar(u), u)
        Df! = (J, du, u, t) -> jacobian!(J, (du, u) -> _f(du, u, t), du, u)
        return RightHandSideFunction(f, f!, Df, Df!)
    end
end
@doc (@doc RightHandSideFunction) RHS(args...; kwargs...) = RightHandSideFunction(args...; kwargs...)