@doc raw"""
    LinearRightHandSide <: AbstractRightHandSide

A composite type for the right-hand side of an [`InitialValueProblem`](@ref) in the form $f(u, t) = L u + g(t)$.

# Constructors
```julia
LinearRightHandSide(L, g, g!)
LinearRightHandSide(L[, g!_or_g]; inplace=nothing)
LRHS(args...; kwargs...)
```

# Arguments
- `L::AbstractMatrix` : $L$, the coefficient term
- `g::Function` : $g$, the forcing term, independent of $u$
- `g!::Function` : $g$ but in-place
"""
struct LinearRightHandSide{
    L_T<:AbstractMatrix{<:Number},g_T<:Union{Function,Nothing},g!_T<:Union{Function,Nothing}
} <: AbstractRightHandSide
    L::L_T
    g::g_T
    g!::g!_T
end

function LinearRightHandSide(L::AbstractMatrix{<:Number}, g!_or_g::Function; inplace=nothing)
    if rhs_inplace(g!_or_g, 1; inplace)
        return LinearRightHandSide(L, _derive_outofplace_forcing(L, g!_or_g), g!_or_g)
    end
    return LinearRightHandSide(L, g!_or_g, _derive_inplace_forcing(g!_or_g))
end
LinearRightHandSide(L::AbstractMatrix{<:Number}) = LinearRightHandSide(L, nothing, nothing)

# A scalar operator becomes a 1×1 matrix.
function LinearRightHandSide(L::Number, g!_or_g::Function; kwargs...)
    return LinearRightHandSide(hcat(L), g!_or_g; kwargs...)
end
LinearRightHandSide(L::Number) = LinearRightHandSide(hcat(L))

@doc (@doc LinearRightHandSide) LRHS(args...; kwargs...) = LinearRightHandSide(args...; kwargs...)

# As for the nonlinear RHS, each derived closure is built in a function of its
# own so that no captured variable is reassigned (and boxed).

# A dense buffer of the operator's element type, whatever the storage of `L`.
function _derive_outofplace_forcing(L, g!)
    return function (t)
        du = Vector{eltype(L)}(undef, size(L, 1))
        g!(du, t)
        return du
    end
end
_derive_inplace_forcing(g) = (du, t) -> du .= g(t)

#----------------------------------- METHODS -----------------------------------

"""
    (rhs::LinearRightHandSide)(u, t)
    (rhs::LinearRightHandSide)(du, u, t)
    (rhs::LinearRightHandSide)(du, v, u, t)

returns the derivative `du` from the solution `u` and time `t`. The 3-argument
in-place form allocates a temporary when a forcing term `g` is present.
For a derived `g(t)`, the dense output buffer has the operator's element type.
Use a complex state in `rhs(u, t)` when a real operator has complex forcing.
The 4-argument form takes caller-owned scratch `v` (same shape as `du`) instead and
is allocation-free. Hot loops should use the 4-argument form.
"""
function (rhs::LinearRightHandSide)(u, t)
    @↓ L, g! = rhs
    du = L * u
    if g! !== nothing
        v = similar(du)
        g!(v, t)
        du .+= v
    end
    return du
end

function (rhs::LinearRightHandSide)(du, u, t)
    @↓ L, g! = rhs
    mul!(du, L, u)
    if g! !== nothing
        v = similar(du)
        g!(v, t)
        du .+= v
    end
    return du
end

function (rhs::LinearRightHandSide)(du, v, u, t)
    @↓ L, g! = rhs
    mul!(du, L, u)
    if g! !== nothing
        g!(v, t)
        du .+= v
    end
    return du
end
