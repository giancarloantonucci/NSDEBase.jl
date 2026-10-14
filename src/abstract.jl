"An abstract type for all objects in NSDEBase."
abstract type AbstractObject end

"""
    SolverFailure <: Exception

An abstract type for numerical failures: a solver could not produce a usable
step (for example, Newton or adaptive-step failure). Custom solvers may subtype
it to request numerical failure handling in Parareal and MoWi. Programming and
configuration errors must use their ordinary exception types instead.
"""
abstract type SolverFailure <: Exception end

"An abstract type for differential-equation problems."
abstract type AbstractProblem <: AbstractObject end

"An abstract type for numerical solvers of [`AbstractProblem`](@ref)s."
abstract type AbstractSolver <: AbstractObject end

"An abstract type for solutions of [`AbstractProblem`](@ref)s."
abstract type AbstractSolution <: AbstractObject end

"An abstract type for parameters used in [`AbstractSolver`](@ref)s."
abstract type AbstractParameters <: AbstractObject end

"An abstract type for caching intermediate computations in [`AbstractSolver`](@ref)s."
abstract type AbstractCache <: AbstractObject end

"""
An abstract type for initial-value problems.

A custom subtype used for chunked or moving-window solves must expose mutable
`rhs`, `u0` and `tspan` fields and implement `remake(problem; u0, tspan)`. The
remade problem must keep the supplied floating-point state buffer: chunk updates
write into it. Solvers also need compatible `initialize_cache`,
`initialize_solution` and `solve!` methods for the custom problem type.
"""
abstract type AbstractInitialValueProblem <: AbstractProblem end

"An abstract type for numerical solvers of [`AbstractInitialValueProblem`](@ref)s."
abstract type AbstractInitialValueSolver <: AbstractSolver end

"An abstract type for computed solutions of [`AbstractInitialValueProblem`](@ref)s."
abstract type AbstractInitialValueSolution <: AbstractSolution end

"An abstract type for parameters in [`AbstractInitialValueSolver`](@ref)s."
abstract type AbstractInitialValueParameters <: AbstractParameters end

"An abstract type for caching intermediate computations in [`AbstractInitialValueSolver`](@ref)s."
abstract type AbstractInitialValueCache <: AbstractCache end

"An abstract type for right-hand-side (RHS) functions of [`AbstractInitialValueProblem`](@ref)s."
abstract type AbstractRightHandSide <: AbstractInitialValueParameters end
