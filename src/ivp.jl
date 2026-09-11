# NSDEBase/src/ivp.jl

"""
    InitialValueProblem <: AbstractInitialValueProblem

A composite type for an initial-value problem.

# Constructors
```julia
InitialValueProblem(rhs, u0, tspan)
InitialValueProblem(rhs, u0, t0, tN)
InitialValueProblem(fₛ, fₙₛ, u0, tspan)
InitialValueProblem(fₛ, fₙₛ, u0, t0, tN)
IVP(args...; kwargs...)
```

# Arguments
- `rhs::AbstractRightHandSide` : right-hand side function
- `u0::AbstractVector` : initial condition
- `tspan::Tuple` : end limits of time domain

# Functions
- [`copy`](@ref) : copy instance of problem with updated parameters
"""
mutable struct InitialValueProblem{rhs_T<:AbstractRightHandSide, u0_T<:AbstractVector{<:Number}, tspan_T<:(Tuple{ℝ,ℝ} where ℝ<:Real)} <: AbstractInitialValueProblem
    rhs::rhs_T
    u0::u0_T
    tspan::tspan_T
end

InitialValueProblem(L::Union{ℂ,AbstractMatrix{ℂ}}, u0::AbstractVector{ℂ}, tspan::Tuple{ℝ,ℝ}) where {ℂ<:Number, ℝ<:Real} = InitialValueProblem(LRHS(L), u0, tspan)
InitialValueProblem(f::Function, u0::AbstractVector{<:Number}, tspan::Tuple{ℝ,ℝ}) where ℝ<:Real = InitialValueProblem(RHS(f), u0, tspan)
InitialValueProblem(rhs::Union{ℂ,AbstractMatrix{ℂ},Function,AbstractRightHandSide}, u0::ℂ, tspan::Tuple{ℝ,ℝ}) where {ℂ<:Number, ℝ<:Real} = InitialValueProblem(rhs, [u0], tspan)
InitialValueProblem(rhs::Union{ℂ,AbstractMatrix{ℂ},Function,AbstractRightHandSide}, u0::Union{ℂ,AbstractVector{ℂ}}, t0::ℝ, tN::ℝ) where {ℂ<:Number, ℝ<:Real} = InitialValueProblem(rhs, u0, (t0, tN))

# SplitRightHandSide construction path, e.g. `IVP(L, fₙₛ, u0, tspan)` for f = L·u + fₙₛ(u, t):
const _SplitStiffPart = Union{Number, AbstractMatrix{<:Number}, Function, LinearRightHandSide, NonlinearRightHandSide}
const _SplitNonStiffPart = Union{Function, NonlinearRightHandSide}
InitialValueProblem(fₛ::_SplitStiffPart, fₙₛ::_SplitNonStiffPart, u0::AbstractVector{<:Number}, tspan::Tuple{ℝ,ℝ}) where ℝ<:Real = InitialValueProblem(SplitRightHandSide(fₛ, fₙₛ), u0, tspan)
InitialValueProblem(fₛ::_SplitStiffPart, fₙₛ::_SplitNonStiffPart, u0::Number, tspan::Tuple{ℝ,ℝ}) where ℝ<:Real = InitialValueProblem(fₛ, fₙₛ, [u0], tspan)
InitialValueProblem(fₛ::_SplitStiffPart, fₙₛ::_SplitNonStiffPart, u0::Union{ℂ,AbstractVector{ℂ}}, t0::ℝ, tN::ℝ) where {ℂ<:Number, ℝ<:Real} = InitialValueProblem(fₛ, fₙₛ, u0, (t0, tN))

@doc (@doc InitialValueProblem) IVP(args...; kwargs...) = InitialValueProblem(args...; kwargs...)

#---------------------------------- FUNCTIONS ----------------------------------

"""
    copy(problem, u0, tspan)::InitialValueProblem
    copy(problem, u0, t0, tN)::InitialValueProblem

returns a copy of `problem` with the same `rhs` but different `u0` and `tspan`.
"""
Base.copy(problem::InitialValueProblem, u0::Union{ℂ,AbstractVector{ℂ}}, tspan::Tuple{ℝ,ℝ}) where {ℂ<:Number, ℝ<:Real} = IVP(problem.rhs, u0, tspan)
Base.copy(problem::InitialValueProblem, u0::Union{ℂ,AbstractVector{ℂ}}, t0::ℝ, tN::ℝ) where {ℂ<:Number, ℝ<:Real} = copy(problem, u0, (t0, tN))
Base.copy(problem::InitialValueProblem, tspan::Tuple{ℝ,ℝ}) where ℝ<:Real = copy(problem, problem.u0, tspan)
Base.copy(problem::InitialValueProblem, t0::ℝ, tN::ℝ) where ℝ<:Real = copy(problem, (t0, tN))
