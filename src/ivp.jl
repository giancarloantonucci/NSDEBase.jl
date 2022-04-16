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
- `rhs :: Union{AbstractRightHandSide, Function, AbstractMatrix{ℂ}, ℂ} where ℂ<:Number`
- `u0 :: Union{AbstractVector{ℂ}, ℂ} where ℂ<:Number`
- `tspan :: Tuple{ℝ, ℝ} where ℝ<:Real`

# Functions
- [`subproblemof`](@ref) : creates a subproblem.
"""
struct InitialValueProblem{rhs_T<:AbstractRightHandSide, u0_T<:(AbstractVector{ℂ} where ℂ<:Number), tspan_T<:(Tuple{ℝ, ℝ} where ℝ<:Real)} <: AbstractInitialValueProblem
    rhs::rhs_T
    u0::u0_T
    tspan::tspan_T
end
InitialValueProblem(L::Union{AbstractMatrix{ℂ}, ℂ}, u0::AbstractVector{ℂ}, tspan::Tuple{ℝ, ℝ}) where {ℂ<:Number, ℝ<:Real} = InitialValueProblem(LRHS(L), u0, tspan)
InitialValueProblem(f::Function, u0::AbstractVector{ℂ}, tspan::Tuple{ℝ, ℝ}) where {ℂ<:Number, ℝ<:Real} = InitialValueProblem(RHS(f), u0, tspan)
InitialValueProblem(rhs::Union{AbstractRightHandSide, Function, AbstractMatrix{ℂ}, ℂ}, u0::ℂ, tspan::Tuple{ℝ, ℝ}) where {ℂ<:Number, ℝ<:Real} = InitialValueProblem(rhs, [u0], tspan)
InitialValueProblem(rhs::Union{AbstractRightHandSide, Function, AbstractMatrix{ℂ}, ℂ}, u0::Union{AbstractVector{ℂ}, ℂ}, t0::ℝ, tN::ℝ) where {ℂ<:Number, ℝ<:Real} = InitialValueProblem(rhs, u0, (t0, tN))
@doc (@doc InitialValueProblem) IVP(args...; kwargs...) = InitialValueProblem(args...; kwargs...)

#####
##### Functions
#####

"""
    subproblemof(problem, u0, tspan) :: InitialValueProblem
    subproblemof(problem, u0, t0, tN) :: InitialValueProblem

returns a subproblem of `problem`, i.e. a copy with same `rhs` but different `u0` and `tspan`.
"""
function subproblemof(problem::InitialValueProblem, u0::Union{ℂ, AbstractVector{ℂ}}, tspan::Tuple{ℝ, ℝ}) where {ℂ<:Number, ℝ<:Real}
    @↓ rhs = problem
    return IVP(rhs, u0, tspan)
end
subproblemof(problem::InitialValueProblem, u0::Union{ℂ, AbstractVector{ℂ}}, t0::ℝ, tN::ℝ) where {ℂ<:Number, ℝ<:Real} = subproblemof(problem, u0, (t0, tN))
