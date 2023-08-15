@doc raw"""
    SplitRightHandSide <: AbstractRightHandSide

A composite type for the right-hand side of an [`InitialValueProblem`](@ref) in the form $f_\text{s}(u, t) + f_\text{ns}(u, t)$.

# Constructors
```julia
SplitRightHandSide(stiff, nonstiff)
SRHS(args...; kwargs...)
```

# Arguments
- `stiff :: Union{LinearRightHandSide, NonlinearRightHandSide}` : $f_\text{s}$, the stiff part of the right-hand side function.
- `nonstiff :: NonlinearRightHandSide}` : $f_\text{ns}$, the non-stiff part of the right-hand side function.
"""
struct SplitRightHandSide{
    stiff_T <: Union{LinearRightHandSide, NonlinearRightHandSide},
    nonstiff_T <: NonlinearRightHandSide,
    } <: AbstractRightHandSide
    stiff :: stiff_T
    nonstiff :: nonstiff_T
end

SplitRightHandSide(f_s::Function, nonstiff::NonlinearRightHandSide) = SplitRightHandSide(NRHS(f_s), nonstiff)
SplitRightHandSide(L::Union{Number, AbstractMatrix{<:Number}}, nonstiff::NonlinearRightHandSide) = SplitRightHandSide(LRHS(L), nonstiff)
SplitRightHandSide(stiff::Union{Number, AbstractMatrix{<:Number}, Function, LinearRightHandSide, NonlinearRightHandSide}, f_ns::Function) = SplitRightHandSide(stiff, NRHS(f_ns))
@doc (@doc SplitRightHandSide) SRHS(args...; kwargs...) = SplitRightHandSide(args...; kwargs...)

#----------------------------------- METHODS -----------------------------------

"""
    (rhs::SplitRightHandSide)(u, t)
    (rhs::SplitRightHandSide)(du, u, t)

returns the derivative `du` from the solution `u` and time `t`.
"""
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
