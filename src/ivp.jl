"""
    InitialValueProblem <: AbstractInitialValueProblem

A composite type for an [`AbstractInitialValueProblem`](@ref).

# Constructors
```julia
InitialValueProblem(rhs, u0, tspan)
InitialValueProblem(rhs, u0, t0, tN)
IVP(args...; kwargs...)
```

# Arguments
- `rhs::Union{Function, AbstractRightHandSideFunction}` : right-hand side derivative.
- `u0::Union{Number, AbstractVector{Number}}` : initial condition.
- `tspan::Tuple{Real, Real}` : time domain.
- `t0::Real` : initial time.
- `tN::Real` : final time.

# Functions
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
struct InitialValueProblem{rhs_T, u0_T, tspan_T} <: AbstractInitialValueProblem
    rhs::rhs_T
    u0::u0_T
    tspan::tspan_T
end

function InitialValueProblem(rhs::AbstractRightHandSideFunction, u0::Number, tspan)
    return InitialValueProblem(rhs, [u0], tspan)
end

function InitialValueProblem(f::Function, u0, tspan)
    return InitialValueProblem(RHS(f), u0, tspan)
end

function InitialValueProblem(rhs::AbstractRightHandSideFunction, u0, t0::Real, tN::Real)
    return InitialValueProblem(rhs, u0, (t0, tN))
end

function InitialValueProblem(f::Function, u0, t0, tN)
    return InitialValueProblem(RHS(f), u0, t0, tN)
end

@doc (@doc InitialValueProblem) function IVP(args...; kwargs...)
    return InitialValueProblem(args...; kwargs...)
end
