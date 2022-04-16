@doc raw"""
    LinearRightHandSide <: AbstractRightHandSide

A composite type for the right-hand side of an [`InitialValueProblem`](@ref) in the form $Lu + g(t)$.

# Constructors
```julia
LinearRightHandSide(L, g, g!)
LinearRightHandSide(L[, g!_or_g])
LRHS(args...; kwargs...)
```

# Arguments
- `L :: Union{Number, AbstractMatrix{<:Number}}`
- `g :: Union{Function, Nothing}`
- `g! :: Union{Function, Nothing}`
"""
struct LinearRightHandSide{L_T<:(AbstractMatrix{ℂ} where ℂ<:Number), g_T<:Union{Function, Nothing}, g!_T<:Union{Function, Nothing}} <: AbstractRightHandSide
    L::L_T
    g::g_T
    g!::g!_T
end
function LinearRightHandSide(L::AbstractMatrix{ℂ}, g!_or_g::Function) where ℂ<:Number
    # Check if g!_or_g is g!(du, t)
    if hasmethod(g!_or_g, NTuple{2, Any})
        g! = g!_or_g
        g = (t) -> g!(similar(L, size(L, 1)), t)
        return LinearRightHandSide(L, g, g!)
    # Check if g!_or_g is g(t)
    elseif hasmethod(g!_or_g, NTuple{1, Any})
        g = g!_or_g
        g! = (du, t) -> du .= g(t)
        return LinearRightHandSide(L, g, g!)
    end
end
LinearRightHandSide(L::Number, g!_or_g::Function) = LinearRightHandSide(hcat(L), g!_or_g)
LinearRightHandSide(L::AbstractMatrix{ℂ}) where ℂ<:Number = LinearRightHandSide(L, nothing, nothing)
LinearRightHandSide(L::Number) = LinearRightHandSide(hcat(L))
@doc (@doc LinearRightHandSide) LRHS(args...; kwargs...) = LinearRightHandSide(args...; kwargs...)

#####
##### Methods
#####

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
    # du = L * u
    mul!(du, L, u)
    if !(g! isa Nothing)
        # du .+= g(t)
        v = similar(du)
        g!(v, t)
        du .+= v
    end
    return du
end
