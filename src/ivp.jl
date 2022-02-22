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
- [`copy`](@ref) : returns a copy.
- [`makesub`](@ref) : creates a subproblem.
"""
struct InitialValueProblem{rhs_T, u0_T, tspan_T} <: AbstractInitialValueProblem
    rhs::rhs_T
    u0::u0_T
    tspan::tspan_T
end

InitialValueProblem(f::Function, u0::AbstractVector, tspan) = InitialValueProblem(RHS(f), u0, tspan)
InitialValueProblem(rhs, u0::Number, tspan) = InitialValueProblem(rhs, [u0], tspan)
InitialValueProblem(rhs, u0, t0, tN) = InitialValueProblem(rhs, u0, (t0, tN))
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
    makesub(problem::InitialValueProblem, u0, tspan) :: InitialValueProblem
    makesub(problem, u0, t0, tN) :: InitialValueProblem

returns a subproblem of `problem`, that is a copy with different `u0` and `tspan`.
"""
function makesub(problem::InitialValueProblem, u0, tspan)
    @↓ rhs = problem
    return IVP(rhs, u0, tspan)
end
makesub(problem, u0, t0, tN) = makesub(problem, u0, (t0, tN))
