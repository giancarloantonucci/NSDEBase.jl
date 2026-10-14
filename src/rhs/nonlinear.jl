@doc raw"""
    NonlinearRightHandSide <: AbstractRightHandSide

A composite type for the right-hand side of an [`InitialValueProblem`](@ref) in the generic form $f(u, t)$.

# Constructors
```julia
NonlinearRightHandSide(f, f!, Df, Df!)
NonlinearRightHandSide(f!_or_f; iscomplex=false, inplace=nothing)
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
`iscomplex = true`. Function-based IVP constructors infer this from complex initial
data. This computes ∂f/∂u; implicit complex solves require a holomorphic RHS.
For non-holomorphic equations, split real and imaginary parts explicitly.
Optional arguments are interpreted by the shortest supported signature; use
`inplace=true` or `inplace=false` to resolve ambiguity.
"""
struct NonlinearRightHandSide{f_T<:Function,f!_T<:Function,Df_T<:Function,Df!_T<:Function} <:
       AbstractRightHandSide
    f::f_T
    f!::f!_T
    Df::Df_T
    Df!::Df!_T
end

# FiniteDifferences uses interleaved real/imaginary coordinates. Extract ∂f/∂u
# directly in O(n²); no dense change-of-basis matrix or linear solve is needed.
function _wirtinger_jacobian(f::Function, u, t)
    J = FiniteDifferences.jacobian(central_fdm(4, 1), u -> f(u, t), u)[1]
    Jrx = @view J[1:2:end, 1:2:end]
    Jry = @view J[1:2:end, 2:2:end]
    Jix = @view J[2:2:end, 1:2:end]
    Jiy = @view J[2:2:end, 2:2:end]
    return complex.((Jrx .+ Jiy) ./ 2, (Jix .- Jry) ./ 2)
end

function rhs_inplace(f, nargs; inplace=nothing)
    inplace === nothing ||
        inplace isa Bool ||
        throw(ArgumentError("inplace must be true, false or nothing"))
    ms = methods(f)
    supports(n) = any(m -> m.isva ? n >= m.nargs - 2 : n == m.nargs - 1, ms)
    signatures = nargs == 2 ? "f(u, t) or f!(du, u, t)" : "g(t) or g!(du, t)"
    if inplace !== nothing
        n = nargs + Int(inplace)
        supports(n) || throw(
            ArgumentError("inplace=$inplace requires a $n-argument method; expected $signatures"),
        )
        return inplace
    end
    arities = [m.nargs - 1 for m in ms if !m.isva]
    nargs in arities && return false
    nargs + 1 in arities && return true
    hint =
        any(m -> m.isva, ms) ? "; choose inplace explicitly for a compatible variadic function" : ""
    return throw(ArgumentError("RHS requires $signatures$hint"))
end

function NonlinearRightHandSide(f!_or_f::Function; iscomplex::Bool=false, inplace=nothing)
    if rhs_inplace(f!_or_f, 2; inplace)
        f!, f = f!_or_f, _derive_outofplace(f!_or_f)
    else
        f!, f = _derive_inplace(f!_or_f), f!_or_f
    end
    Df, Df! = iscomplex ? _wirtinger_jacobians(f) : _forwarddiff_jacobians(f, f!)
    return NonlinearRightHandSide(f, f!, Df, Df!)
end

@doc (@doc NonlinearRightHandSide) RightHandSide(args...; kwargs...) =
    NonlinearRightHandSide(args...; kwargs...)
@doc (@doc RightHandSide) RHS(args...; kwargs...) = RightHandSide(args...; kwargs...)

# Each derived closure is built in a function of its own, so no captured
# variable is ever reassigned: Julia boxes such variables, which makes every
# call through the closure type-unstable.

# Return the filled buffer, never `f!`'s own return value: a mutating function
# is entitled to return `nothing`, and the out-of-place path (and ForwardDiff
# through it) needs the array.
_derive_outofplace(f!) = (u, t) -> (du=similar(u); f!(du, u, t); du)
_derive_inplace(f) = (du, u, t) -> du .= f(u, t)

function _forwarddiff_jacobians(f, f!)
    Df = (u, t) -> ForwardDiff.jacobian(u -> f(u, t), u)
    Df! = (J, du, u, t) -> ForwardDiff.jacobian!(J, (du, u) -> f!(du, u, t), du, u)
    return Df, Df!
end

function _wirtinger_jacobians(f)
    Df = (u, t) -> _wirtinger_jacobian(f, u, t)
    Df! = (J, du, u, t) -> J .= Df(u, t)
    return Df, Df!
end

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
    f!(du, u, t)
    return du
end

(rhs::NonlinearRightHandSide)(du, v, u, t) = rhs(du, u, t)
