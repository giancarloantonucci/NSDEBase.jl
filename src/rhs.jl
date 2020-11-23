@doc """
    RightHandSideFunction(f[, Df]) -> RightHandSideFunction
    RHS(args...; kwargs...) -> RightHandSideFunction

returns a constructor for the right-hand side derivative of an `InitialValueProblem` with
- `f` : right-hand side derivative (standard or in-place).
- `Df` : jacobian of `f` (standard or in-place).
"""
struct RightHandSideFunction{f_T, Df_T}
    f::f_T
    Df::Df_T
end

function RightHandSideFunction(_f::Function)
    if hasmethod(_f, NTuple{2,Any})
        f = (u, t) -> _f(u, t)
        Df = (u, t) -> jacobian(u -> f(u, t), u)
        return RightHandSideFunction(f, Df)
    elseif hasmethod(_f, NTuple{3,Any})
        f = (du, u, t) -> _f(du, u, t)
        Df = (J, u, t) -> jacobian!(J, (du, u) -> f(du, u, t), similar(u), u)
        return RightHandSideFunction(f, Df)
    end
end
@doc (@doc RightHandSideFunction) RHS(args...; kwargs...) = RightHandSideFunction(args...; kwargs...)
