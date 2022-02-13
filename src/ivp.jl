"""
    InitialValueProblem <: AbstractInitialValueProblem

A composite type for an [`AbstractInitialValueProblem`](@ref).

# Constructors
```julia
InitialValueProblem(rhs, u0, tspan)
InitialValueProblem(rhs, u0, t0, tN)
IVP(args...; kwargs...)
```

## Arguments
- `rhs :: Union{Function, AbstractRightHandSideFunction}` : right-hand side derivative.
- `u0 :: Union{Number, AbstractVector{<:Number}}` : initial condition.
- `tspan :: Tuple{Real, Real}` : time domain.
- `t0 :: Real` : initial time.
- `tN :: Real` : final time.

# Functions
- [`subproblemof`](@ref) : creates a subproblem.
"""
struct InitialValueProblem{rhs_T, u0_T, tspan_T} <: AbstractInitialValueProblem
    rhs::rhs_T
    u0::u0_T
    tspan::tspan_T
end

InitialValueProblem(rhs::AbstractRightHandSideFunction, u0::AbstractVector{<:AbstractFloat}, tspan::Tuple{Integer, Integer}) = InitialValueProblem(rhs, u0, float.(tspan))
InitialValueProblem(rhs::AbstractRightHandSideFunction, u0::AbstractVector{<:Integer}, tspan) = InitialValueProblem(rhs, float.(u0), tspan)
InitialValueProblem(rhs::AbstractRightHandSideFunction, u0::Number, tspan) = InitialValueProblem(rhs, [u0], tspan)
InitialValueProblem(f::Function, u0, tspan) = InitialValueProblem(RHS(f), u0, tspan)
InitialValueProblem(rhs::AbstractRightHandSideFunction, u0, t0, tN) = InitialValueProblem(rhs, u0, (t0, tN))
InitialValueProblem(f::Function, u0, t0, tN) = InitialValueProblem(RHS(f), u0, t0, tN)
@doc (@doc InitialValueProblem) IVP(args...; kwargs...) = InitialValueProblem(args...; kwargs...)

#####
##### Functions
#####

"""
    copy(problem::InitialValueProblem)
    
returns a copy of `problem`.
"""
function Base.copy(problem::InitialValueProblem)
    @↓ rhs, u0, tspan = problem
    return IVP(rhs, u0, tspan)
end

"""
    subproblemof(problem::InitialValueProblem, u0, tspan) :: InitialValueProblem
    subproblemof(problem, u0, t0, tN) :: InitialValueProblem

returns a copy of `problem` with different initial condition and time domain.
"""
function subproblemof(problem::InitialValueProblem, u0, tspan)
    @↓ rhs = problem
    return IVP(rhs, u0, tspan)
end
subproblemof(problem, u0, t0, tN) = subproblemof(problem, u0, (t0, tN))
