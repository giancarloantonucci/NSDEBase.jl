# NSDEBase/src/ivp.jl

"""
    InitialValueProblem <: AbstractInitialValueProblem

A composite type for an initial-value problem. Time endpoints are promoted to a
common floating-point type; integer and rational states are converted to floating
point. Already-floating state vectors (including complex states) retain their
storage. The solver time grid uses the resulting problem time type.

# Constructors
```julia
InitialValueProblem(rhs, u0, tspan)
InitialValueProblem(rhs, u0, t0, tN)
InitialValueProblem(fₛ, fₙₛ, u0, tspan)
InitialValueProblem(fₛ, fₙₛ, u0, t0, tN)
IVP(args...; kwargs...)
```

Function constructors accept `inplace=true/false/nothing` for scalar or vector
states and tuple or separate time endpoints. Split constructors also accept
`inplace=(stiff_choice, nonstiff_choice)`; a shared choice applies to both raw
functions. Prebuilt RHS objects retain their own conventions.

# Arguments
- `rhs::AbstractRightHandSide` : right-hand side function
- `u0::AbstractVector` : initial condition
- `tspan::Tuple` : end limits of time domain

# Functions
- [`remake`](@ref): construct a problem with selected updated parameters
"""
mutable struct InitialValueProblem{
    rhs_T<:AbstractRightHandSide,u0_T<:AbstractVector{<:Number},tspan_T<:(Tuple{ℝ,ℝ} where ℝ<:Real)
} <: AbstractInitialValueProblem
    rhs::rhs_T
    u0::u0_T
    tspan::tspan_T
    function InitialValueProblem(
        rhs::AbstractRightHandSide, u0::AbstractVector{<:Number}, tspan::Tuple{Real,Real}
    )
        ts = float.(promote(tspan...))
        T = state_eltype(rhs, float(eltype(u0)))
        v = eltype(u0) === T ? u0 : T.(u0)
        return new{typeof(rhs),typeof(v),typeof(ts)}(rhs, v, ts)
    end
end

state_eltype(rhs::AbstractRightHandSide, T) = T
function state_eltype(rhs::LinearRightHandSide, T)
    return eltype(rhs.L) <: Complex ? promote_type(T, eltype(rhs.L)) : T
end
state_eltype(rhs::SplitRightHandSide, T) = state_eltype(rhs.fₛ, T)

function InitialValueProblem(
    L::Union{Number,AbstractMatrix{<:Number}}, u0::AbstractVector{<:Number}, tspan::Tuple{Real,Real}
)
    return InitialValueProblem(LRHS(L), u0, tspan)
end
function InitialValueProblem(
    f::Function, u0::AbstractVector{<:Number}, tspan::Tuple{Real,Real}; inplace=nothing
)
    return InitialValueProblem(RHS(f; iscomplex=eltype(u0) <: Complex, inplace), u0, tspan)
end
function InitialValueProblem(
    rhs::Union{Number,AbstractMatrix{<:Number},Function,AbstractRightHandSide},
    u0::Number,
    tspan::Tuple{Real,Real};
    kwargs...,
)
    return InitialValueProblem(rhs, [u0], tspan; kwargs...)
end
function InitialValueProblem(
    rhs::Union{Number,AbstractMatrix{<:Number},Function,AbstractRightHandSide},
    u0::Union{Number,AbstractVector{<:Number}},
    t0::Real,
    tN::Real;
    kwargs...,
)
    return InitialValueProblem(rhs, u0, (t0, tN); kwargs...)
end

# SplitRightHandSide construction path, e.g. `IVP(L, fₙₛ, u0, tspan)` for f = L·u + fₙₛ(u, t):
const _SplitStiffPart = Union{
    Number,AbstractMatrix{<:Number},Function,LinearRightHandSide,NonlinearRightHandSide
}
const _SplitNonStiffPart = Union{Function,NonlinearRightHandSide}
function InitialValueProblem(
    fₛ::_SplitStiffPart,
    fₙₛ::_SplitNonStiffPart,
    u0::AbstractVector{<:Number},
    tspan::Tuple{Real,Real};
    inplace=nothing,
)
    stiff_inplace, nonstiff_inplace = split_inplace(inplace)
    stiff = fₛ isa Union{Number,AbstractMatrix} ? LRHS(fₛ) : fₛ
    T = stiff isa AbstractRightHandSide ? state_eltype(stiff, float(eltype(u0))) : float(eltype(u0))
    stiff = stiff isa Function ? RHS(stiff; iscomplex=T <: Complex, inplace=stiff_inplace) : stiff
    nonstiff = fₙₛ isa Function ? RHS(fₙₛ; iscomplex=T <: Complex, inplace=nonstiff_inplace) : fₙₛ
    return InitialValueProblem(SplitRightHandSide(stiff, nonstiff), u0, tspan)
end
function InitialValueProblem(
    fₛ::_SplitStiffPart, fₙₛ::_SplitNonStiffPart, u0::Number, tspan::Tuple{Real,Real}; kwargs...
)
    return InitialValueProblem(fₛ, fₙₛ, [u0], tspan; kwargs...)
end
function InitialValueProblem(
    fₛ::_SplitStiffPart,
    fₙₛ::_SplitNonStiffPart,
    u0::Union{Number,AbstractVector{<:Number}},
    t0::Real,
    tN::Real;
    kwargs...,
)
    return InitialValueProblem(fₛ, fₙₛ, u0, (t0, tN); kwargs...)
end

@doc (@doc InitialValueProblem) IVP(args...; kwargs...) = InitialValueProblem(args...; kwargs...)

#---------------------------------- FUNCTIONS ----------------------------------

"""
    remake(problem; u0=problem.u0, tspan=problem.tspan)

Build a problem with the same right-hand side and selected new initial data.
Time-span endpoints are promoted to a common numeric type.
A supplied floating-point state is not copied; other numeric states are converted.
Use `copy(u0)` for independent storage.
"""
function remake(problem::InitialValueProblem; u0=problem.u0, tspan=problem.tspan)
    return IVP(problem.rhs, u0, promote(tspan...))
end

"""
    chunk_solver(solver, duration)

Prepare a solver for a chunk of nominal duration. Solver packages may specialise
this interface to share a fixed schedule across translated chunks. The fallback
copies the solver, so run-local mutable state is never shared with its caller.
"""
chunk_solver(solver::AbstractInitialValueSolver, duration::Real) = deepcopy(solver)

"Return a solver with changed immutable step-size settings; implemented by solver packages."
function withstepsize end

function remake(problem::AbstractInitialValueProblem; kwargs...)
    return throw(
        ArgumentError(
            "$(typeof(problem)) must implement remake(problem; u0, tspan) for chunked and moving-window solves",
        ),
    )
end
