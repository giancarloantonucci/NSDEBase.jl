"""
    SplitRightHandSideFunction <: AbstractRightHandSideFunction

A composite type for the split right-hand side of an [`AbstractInitialValueProblem`](@ref).

# Constructors
```julia
SplitRightHandSideFunction(L, rhs)
SRHS(args...; kwargs...)
```

# Arguments
- `L :: Union{Number, Matrix}` : linear part of right-hand side derivative.
- `rhs :: RightHandSideFunction` : nonlinear part of right-hand side derivative.
"""
mutable struct SplitRightHandSideFunction{L_T, rhs_T} <: AbstractRightHandSideFunction
    L::L_T
    rhs::rhs_T
end

@doc (@doc SplitRightHandSideFunction) SRHS(args...; kwargs...) = SplitRightHandSideFunction(args...; kwargs...)
