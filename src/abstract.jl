"An abstract type for all objects in NSDEBase."
abstract type AbstractObject end

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

"An abstract type for initial value problems (IVPs)."
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
