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
- `rhs   :: Union{Function, AbstractRightHandSideFunction}` : right-hand side derivative.
- `u0    :: Union{Number, AbstractVector{<:Number}}`        : initial condition.
- `tspan :: Tuple{Real, Real}`                              : time domain.
- `t0    :: Real`                                           : initial time.
- `tN    :: Real`                                           : final time.

# Functions
- [`show`](@ref)         : shows name and contents.
- [`subproblemof`](@ref) : creates a subproblem.
- [`summary`](@ref)      : shows name.
"""
struct InitialValueProblem{rhs_T<:Union{Function, AbstractRightHandSideFunction}, u0_T<:Union{Number, AbstractVector{<:Number}}, tspan_T<:Tuple{Real, Real}} <: AbstractInitialValueProblem
    rhs::rhs_T
    u0::u0_T
    tspan::tspan_T
end

InitialValueProblem(rhs::AbstractRightHandSideFunction, u0::Number, tspan) = InitialValueProblem(rhs, [u0], tspan)
InitialValueProblem(f::Function, u0, tspan) = InitialValueProblem(RHS(f), u0, tspan)
InitialValueProblem(rhs::AbstractRightHandSideFunction, u0, t0, tN) = InitialValueProblem(rhs, u0, (t0, tN))
InitialValueProblem(f::Function, u0, t0, tN) = InitialValueProblem(RHS(f), u0, t0, tN)
@doc (@doc InitialValueProblem) IVP(args...; kwargs...) = InitialValueProblem(args...; kwargs...)

#####
##### Functions
#####

"""
    subproblemof(problem::InitialValueProblem, u0, tspan) :: InitialValueProblem
    subproblemof(problem, u0, t0, tN) :: InitialValueProblem

returns a copy of `problem` with different initial condition and time domain.
"""
subproblemof(problem::InitialValueProblem, u0, tspan) = IVP(problem.rhs, u0, tspan)
subproblemof(problem, u0, t0, tN) = subproblemof(problem, u0, (t0, tN))
