@doc raw"""
    SplitRightHandSide <: AbstractRightHandSide

A composite type for the right-hand side of an [`InitialValueProblem`](@ref) in the form $f_s(u, t) + f_{ns}(u, t)$.

# Constructors
```julia
SplitRightHandSide(stiff, nonstiff)
SRHS(args...; kwargs...)
```

# Arguments
- `stiff :: Union{NonlinearRightHandSide, LinearRightHandSide, Function, AbstractMatrix{ℂ}, ℂ} where ℂ<:Number` : $f_s$.
- `nonstiff :: Union{NonlinearRightHandSide, Function}` : $f_{ns}$.
"""
struct SplitRightHandSide{stiff_T<:Union{LinearRightHandSide, NonlinearRightHandSide}, nonstiff_T<:NonlinearRightHandSide} <: AbstractRightHandSide
    stiff::stiff_T
    nonstiff::nonstiff_T
end

SplitRightHandSide(fₛ::Function, nonstiff::NonlinearRightHandSide) = SplitRightHandSide(NRHS(fₛ), nonstiff)
SplitRightHandSide(L::Union{AbstractMatrix{ℂ}, ℂ}, nonstiff::NonlinearRightHandSide) where ℂ<:Number = SplitRightHandSide(LRHS(L), nonstiff)
SplitRightHandSide(stiff::Union{NonlinearRightHandSide, LinearRightHandSide, Function, AbstractMatrix{ℂ}, ℂ}, fₙₛ::Function) where ℂ<:Number = SplitRightHandSide(stiff, NRHS(fₙₛ))
@doc (@doc SplitRightHandSide) SRHS(args...; kwargs...) = SplitRightHandSide(args...; kwargs...)

#----------------------------------- METHODS -----------------------------------

function (rhs::SplitRightHandSide)(u, t)
    @↓ stiff, nonstiff = rhs
    du = stiff(u, t) + nonstiff(u, t)
    return du
end

function (rhs::SplitRightHandSide)(du, u, t)
    @↓ stiff, nonstiff = rhs
    # @! du = stiff(u, t) + nonstiff(u, t)
    stiff(du, u, t)
    v = similar(du)
    nonstiff(v, u, t)
    du .+= v
    return du
end
