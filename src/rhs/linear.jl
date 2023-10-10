@doc raw"""
    LinearRightHandSide <: AbstractRightHandSide

A composite type for the right-hand side of an [`InitialValueProblem`](@ref) in the form $f(u, t) = L u + g(t)$.

# Constructors
```julia
LinearRightHandSide(L, g, g!)
LinearRightHandSide(L[, g!_or_g])
LRHS(args...; kwargs...)
```

# Arguments
- `L::AbstractMatrix` : $L$, the coefficient term
- `g::Function` : $g$, the forcing term, independent of $u$
- `g!::Function` : $g$ (in-place)
"""
struct LinearRightHandSide{L_T<:AbstractMatrix{<:Number}, g_T<:Union{Function,Nothing}, g!_T<:Union{Function,Nothing}} <: AbstractRightHandSide
    L::L_T
    g::g_T
    g!::g!_T
end

function LinearRightHandSide(L::AbstractMatrix{<:Number}, g!_or_g::Function)
    if hasmethod(g!_or_g, NTuple{2, Any}) # has g!_or_g signature g!(du, t)?
        g! = g!_or_g
        g = t -> g!(similar(L, size(L, 1)), t)
        return LinearRightHandSide(L, g, g!)
    elseif hasmethod(g!_or_g, NTuple{1, Any}) # has g!_or_g signature g(t)?
        g = g!_or_g
        g! = (du, t) -> du .= g(t)
        return LinearRightHandSide(L, g, g!)
    else
        throw(ArgumentError("`LinearRightHandSide(L, g!_or_g)` needs `g!_or_g` to have signature `g!(du, t)` or `g(t)`."))
    end
end

LinearRightHandSide(L::Number, g!_or_g::Function) = LinearRightHandSide(hcat(L), g!_or_g) # hcat(::Number)::Matrix
LinearRightHandSide(L::AbstractMatrix{<:Number}) = LinearRightHandSide(L, nothing, nothing)
LinearRightHandSide(L::Number) = LinearRightHandSide(hcat(L)) # hcat(Number)::Matrix
@doc (@doc LinearRightHandSide) LRHS(args...; kwargs...) = LinearRightHandSide(args...; kwargs...)

#----------------------------------- METHODS -----------------------------------

"""
    (rhs::LinearRightHandSide)(u, t)
    (rhs::LinearRightHandSide)(du, u, t)

returns the derivative `du` from the solution `u` and time `t`.
"""
function (rhs::LinearRightHandSide)(u, t)
    @↓ L, g = rhs
    du = L * u
    if !(g isa Nothing)
        du .+= g(t)
    end
    return du
end

function (rhs::LinearRightHandSide)(du, u, t)
    @↓ L, g! = rhs
    # @! du = L * u
    mul!(du, L, u)
    if !(g! isa Nothing)
        # @! du .+= g(t)
        v = similar(du)
        g!(v, t)
        du .+= v
    end
    return du
end
