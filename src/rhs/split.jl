# NSDEBase/src/rhs/split.jl

@doc raw"""
    SplitRightHandSide <: AbstractRightHandSide

A composite type for the right-hand side of an [`InitialValueProblem`](@ref) in the form $f(u, t) = f_\text{s}(u, t) + f_\text{ns}(u, t)$.

# Constructors
```julia
SplitRightHandSide(fₛ, fₙₛ)
SRHS(args...; kwargs...)
```

# Arguments
- `fₛ::Union{LinearRightHandSide,NonlinearRightHandSide}` : $f_\text{s}$, the stiff part of the right-hand side function $f$
- `fₙₛ::NonlinearRightHandSide` : $f_\text{ns}$, the non-stiff part of the right-hand side function $f$
"""
struct SplitRightHandSide{fₛ_T<:Union{LinearRightHandSide,NonlinearRightHandSide}, fₙₛ_T<:NonlinearRightHandSide} <: AbstractRightHandSide
    fₛ::fₛ_T
    fₙₛ::fₙₛ_T
end

SplitRightHandSide(fₛ::Function, fₙₛ::NonlinearRightHandSide) = SplitRightHandSide(NonlinearRightHandSide(fₛ), fₙₛ)
SplitRightHandSide(L::Union{Number,AbstractMatrix{<:Number}}, fₙₛ::NonlinearRightHandSide) = SplitRightHandSide(LinearRightHandSide(L), fₙₛ)
SplitRightHandSide(fₛ::Union{Number,AbstractMatrix{<:Number},Function,LinearRightHandSide,NonlinearRightHandSide}, fₙₛ::Function) = SplitRightHandSide(fₛ, NonlinearRightHandSide(fₙₛ))
@doc (@doc SplitRightHandSide) SRHS(args...; kwargs...) = SplitRightHandSide(args...; kwargs...)

#----------------------------------- METHODS -----------------------------------

"""
    (rhs::SplitRightHandSide)(u, t)
    (rhs::SplitRightHandSide)(du, u, t)
    (rhs::SplitRightHandSide)(du, v, u, t)

returns the derivative `du` from the solution `u` and time `t`. The 3-argument
in-place form allocates a temporary; the 4-argument form takes caller-owned
scratch `v` (same shape as `du`) instead and is allocation-free. `v` is used
first as scratch for `fₛ` and then to hold `fₙₛ(u, t)`, so a single buffer
covers the whole evaluation. Hot loops should use the 4-argument form.
"""
function (rhs::SplitRightHandSide)(u, t)
    @↓ fₛ, fₙₛ = rhs
    du = fₛ(u, t) + fₙₛ(u, t)
    return du
end

function (rhs::SplitRightHandSide)(du, u, t)
    @↓ fₛ, fₙₛ = rhs
    # @! du = fₛ(u, t) + fₙₛ(u, t)
    fₛ(du, u, t)
    v = similar(du)
    fₙₛ(v, u, t)
    du .+= v
    return du
end

function (rhs::SplitRightHandSide)(du, v, u, t)
    @↓ fₛ, fₙₛ = rhs
    # @! du = fₛ(u, t) + fₙₛ(u, t)
    fₛ(du, v, u, t) # fₛ is done with the scratch `v` once it returns
    fₙₛ(v, u, t)
    du .+= v
    return du
end
