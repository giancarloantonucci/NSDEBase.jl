# NSDEBase/src/rhs/nonlinear.jl

@doc raw"""
    NonlinearRightHandSide <: AbstractRightHandSide

A composite type for the right-hand side of an [`InitialValueProblem`](@ref) in the generic form $f(u, t)$.

# Constructors
```julia
NonlinearRightHandSide(f, f!, Df, Df!)
NonlinearRightHandSide(f!_or_f; iscomplex=false)
RightHandSide(args...; kwargs...)
RHS(args...; kwargs...)
```

# Arguments
- `f::Function` : $f$, the right-hand side function
- `f!::Function` : $f$ but in-place
- `Df::Function` : $\mathcal{D}f$, the Jacobian of $f$ with respect to $u$
- `Df!::Function` : $\mathcal{D}f$ but in-place

When only `f!_or_f` is given, the missing form is derived from it and the
Jacobian is filled in automatically: ForwardDiff for real problems, and a
finite-difference Wirtinger derivative (via FiniteDifferences) when
`iscomplex = true`.
"""
struct NonlinearRightHandSide{f_T<:Function, f!_T<:Function, Df_T<:Function, Df!_T<:Function} <: AbstractRightHandSide
    f::f_T
    f!::f!_T
    Df::Df_T
    Df!::Df!_T
end

# Wirtinger derivative ∂f/∂u of a ℂⁿ → ℂⁿ map, via a real finite-difference
# Jacobian in the (Re, Im) basis followed by a change of basis to (u, conj(u)).
function _wirtinger_jacobian(f::Function, u, t)
    J = FiniteDifferences.jacobian(central_fdm(4, 1), u -> f(u, t), u)[1]
    n = length(u)
    I_n = Matrix(I, n, n)
    P = kron(I_n, 0.5 * [1.0 -1.0im; 1.0 1.0im])
    J_p = (P * J) / P # switch basis
    # Per element, Q = ½[1 -i; 1 i] maps [Re(f), Im(f)] ↦ ½[conj(f), f]: row 1
    # of each 2×2 block is the conjugate output, row 2 is f itself (likewise
    # for the input columns). So the ∂f/∂u block sits at the EVEN indices;
    # the odd-index block is ∂conj(f)/∂conj(u) = conj(∂f/∂u).
    return J_p[2:2:2n, 2:2:2n] # keep the ∂f/∂u block
end

function NonlinearRightHandSide(f!_or_f::Function; iscomplex::Bool=false)
    if hasmethod(f!_or_f, NTuple{3, Any}) # i.e. has f!_or_f signature f!(du, u, t)?
        f! = f!_or_f
        # Return the filled buffer, never `f!`'s own return value: a mutating
        # function is entitled to return `nothing`, and the out-of-place path
        # (and ForwardDiff through it) needs the array.
        f = (u, t) -> (du = similar(u); f!(du, u, t); du)
    elseif hasmethod(f!_or_f, NTuple{2, Any}) # i.e. has f!_or_f signature f(u, t)?
        f = f!_or_f
        f! = (du, u, t) -> du .= f(u, t)
    else
        throw(ArgumentError("`NonlinearRightHandSide(f!_or_f; ...)` needs `f!_or_f` to have signature `f!(du, u, t)` or `f(u, t)`."))
    end
    if iscomplex
        Df = (u, t) -> _wirtinger_jacobian(f, u, t)
        Df! = (J, du, u, t) -> J .= Df(u, t)
    else
        Df = (u, t) -> ForwardDiff.jacobian(u -> f(u, t), u)
        Df! = (J, du, u, t) -> ForwardDiff.jacobian!(J, (du, u) -> f!(du, u, t), du, u)
    end
    return NonlinearRightHandSide(f, f!, Df, Df!)
end

@doc (@doc NonlinearRightHandSide) RightHandSide(args...; kwargs...) = NonlinearRightHandSide(args...; kwargs...)
@doc (@doc RightHandSide) RHS(args...; kwargs...) = RightHandSide(args...; kwargs...)

#----------------------------------- METHODS -----------------------------------

"""
    (rhs::NonlinearRightHandSide)(u, t)
    (rhs::NonlinearRightHandSide)(du, u, t)
    (rhs::NonlinearRightHandSide)(du, v, u, t)

returns the derivative `du` from the solution `u` and time `t`. The 4-argument
form takes caller-owned scratch `v` (same shape as `du`) so that hot loops can
guarantee an allocation-free call; for a nonlinear RHS the scratch is unused.
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

(rhs::NonlinearRightHandSide)(du, v, u, t) = rhs(du, u, t)
