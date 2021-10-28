"""
    SplitRightHandSideFunction <: AbstractRightHandSideFunction

A composite type for the split right-hand side of an [`AbstractInitialValueProblem`](@ref).

# Constructors
```julia
SplitRightHandSideFunction(L, rhs)
SRHS(args...; kwargs...)
```

# Arguments
- `L::Union{Number, AbstractMatrix}` : linear part of right-hand side derivative.
- `rhs::RightHandSideFunction` : nonlinear part of right-hand side derivative.

# Functions
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
struct SplitRightHandSideFunction{L_T, rhs_T} <: AbstractRightHandSideFunction
    L::L_T
    rhs::rhs_T
end

@doc (@doc SplitRightHandSideFunction) function SRHS(args...; kwargs...)
    return SplitRightHandSideFunction(args...; kwargs...)
end
