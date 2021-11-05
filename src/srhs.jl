"""
    SplitRightHandSideFunction <: AbstractRightHandSideFunction

A composite type for the split right-hand side of an [`InitialValueProblem`](@ref).

# Constructors
```julia
SplitRightHandSideFunction(L, rhs)
SRHS(args...; kwargs...)
```

## Arguments
- `L :: Union{Number, AbstractMatrix{<:Number}}` : linear part of right-hand side derivative.
- `rhs :: RightHandSideFunction` : nonlinear part of right-hand side derivative.
"""
struct SplitRightHandSideFunction{L_T<:Union{Number, AbstractMatrix{<:Number}}, rhs_T<:RightHandSideFunction} <: AbstractRightHandSideFunction
    L::L_T
    rhs::rhs_T
end

@doc (@doc SplitRightHandSideFunction) SRHS(args...; kwargs...) = SplitRightHandSideFunction(args...; kwargs...)
