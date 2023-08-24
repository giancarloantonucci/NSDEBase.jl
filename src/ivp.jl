"""
    InitialValueProblem <: AbstractInitialValueProblem

A composite type for an initial-value problem.

# Constructors
```julia
InitialValueProblem(rhs, u0, tspan)
InitialValueProblem(rhs, u0, t0, tN)
IVP(args...; kwargs...)
```

# Arguments
- `rhs :: AbstractRightHandSide` : right-hand side function.
- `u0 :: AbstractVector` : initial condition.
- `tspan :: Tuple` : end limits of time domain.

# Functions
- [`new`](@ref) : new instance of problem with updated parameters.
"""
struct InitialValueProblem{
    rhs_T   <: AbstractRightHandSide,
    u0_T    <: AbstractVector{<:Number},
    tspan_T <: (Tuple{ℝ, ℝ} where ℝ <: Real),
    } <: AbstractInitialValueProblem
    rhs   :: rhs_T
    u0    :: u0_T
    tspan :: tspan_T
end

# TODO: Add methods for SRHS
InitialValueProblem(L::Union{ℂ, AbstractMatrix{ℂ}}, u0::AbstractVector{ℂ}, tspan::Tuple{ℝ, ℝ}) where {ℂ <: Number, ℝ <: Real} = InitialValueProblem(LRHS(L), u0, tspan)
InitialValueProblem(f::Function, u0::AbstractVector{<:Number}, tspan::Tuple{ℝ, ℝ}) where ℝ <: Real = InitialValueProblem(RHS(f), u0, tspan)
InitialValueProblem(rhs::Union{ℂ, AbstractMatrix{ℂ}, Function, AbstractRightHandSide}, u0::ℂ, tspan::Tuple{ℝ, ℝ}) where {ℂ <: Number, ℝ <: Real} = InitialValueProblem(rhs, [u0], tspan)
InitialValueProblem(rhs::Union{ℂ, AbstractMatrix{ℂ}, Function, AbstractRightHandSide}, u0::Union{ℂ, AbstractVector{ℂ}}, t0::ℝ, tN::ℝ) where {ℂ <: Number, ℝ <: Real} = InitialValueProblem(rhs, u0, (t0, tN))
@doc (@doc InitialValueProblem) IVP(args...; kwargs...) = InitialValueProblem(args...; kwargs...)

#---------------------------------- FUNCTIONS ----------------------------------

"""
    new(problem, u0, tspan) :: InitialValueProblem
    new(problem, u0, t0, tN) :: InitialValueProblem

returns a new instance of `problem` with the same `rhs` but different `u0` and `tspan`.
"""
new(problem::InitialValueProblem, u0::Union{ℂ, AbstractVector{ℂ}}, tspan::Tuple{ℝ, ℝ}) where {ℂ <: Number, ℝ <: Real} = IVP(problem.rhs, u0, tspan)
new(problem::InitialValueProblem, u0::Union{ℂ, AbstractVector{ℂ}}, t0::ℝ, tN::ℝ) where {ℂ <: Number, ℝ <: Real} = new(problem, u0, (t0, tN))
