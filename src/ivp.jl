"An abstract type for the solvers of [`InitialValueProblem`](@ref)s."
abstract type InitialValueSolver <: NSDESolver end

"An abstract type for the solution of an [`InitialValueProblem`](@ref)."
abstract type InitialValueSolution <: NSDESolution end

"""
    InitialValueProblem <: NSDEProblem

A composite type for initial value problems.

# Constructors
```julia
InitialValueProblem(rhs, u0, tspan)
InitialValueProblem(rhs, u0, t0, tN)
IVP(args...; kwargs...)
```

# Arguments
- `rhs :: Union{Function, RightHandSideFunction}` : right-hand side derivative.
- `u0 :: Union{Number, AbstractVector{Number}}` : initial condition.
- `tspan :: Tuple{Real, Real}` : time domain.
- `t0 :: Real` : initial time.
- `tN :: Real` : final time.

# Functions
- [`show`](@ref) : shows name and contents.
- [`summary`](@ref) : shows name.
"""
struct InitialValueProblem{rhs_T, u0_T, tspan_T} <: NSDEProblem
    rhs::rhs_T
    u0::u0_T
    tspan::tspan_T
end

InitialValueProblem(rhs::RightHandSideFunction, u0::Number, tspan) = InitialValueProblem(rhs, [u0], tspan)
InitialValueProblem(f::Function, u0, tspan) = InitialValueProblem(RHS(f), u0, tspan)
InitialValueProblem(rhs::RightHandSideFunction, u0, t0::Real, tN::Real) = InitialValueProblem(rhs, u0, (t0, tN))
InitialValueProblem(f::Function, u0, t0, tN) = InitialValueProblem(RHS(f), u0, t0, tN)
@doc (@doc InitialValueProblem) IVP(args...; kwargs...) = InitialValueProblem(args...; kwargs...)

# ---------------------------------------------------------------------------- #
#                                   Functions                                  #
# ---------------------------------------------------------------------------- #

"""
    show(io::IO, problem::InitialValueProblem)

prints a full description of `problem` and its contents to a stream `io`.
"""
Base.show(io::IO, problem::InitialValueProblem) = _show(io, problem)

"""
    summary(io::IO, problem::InitialValueProblem)

prints a brief description of `problem` to a stream `io`.
"""
Base.summary(io::IO, problem::InitialValueProblem) = _summary(io, problem)
