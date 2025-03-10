"""
    AbstractObject

Abstract type for all objects in NSDEBase.
"""
abstract type AbstractObject end

"""
    AbstractProblem <: AbstractObject

Abstract type for differential-equation problems.
"""
abstract type AbstractProblem <: AbstractObject end

"""
    AbstractSolver <: AbstractObject

Abstract type for numerical solvers of [`AbstractProblem`](@ref)s.
"""
abstract type AbstractSolver <: AbstractObject end

"""
    AbstractSolution <: AbstractObject

Abstract type for solutions of [`AbstractProblem`](@ref)s.
"""
abstract type AbstractSolution <: AbstractObject end

"""
    AbstractParameters <: AbstractObject

Abstract type for parameter containers used in [`AbstractSolver`](@ref)s.
"""
abstract type AbstractParameters <: AbstractObject end

"""
    AbstractCache <: AbstractObject

Abstract type for caching intermediate computations in [`AbstractSolver`](@ref)s.
"""
abstract type AbstractCache <: AbstractObject end

"""
    AbstractInitialValueProblem <: AbstractProblem

Abstract type for initial value problems (IVPs).
"""
abstract type AbstractInitialValueProblem <: AbstractProblem end

"""
    AbstractInitialValueSolver <: AbstractSolver

Abstract type for numerical solvers of [`AbstractInitialValueProblem`](@ref)s.
"""
abstract type AbstractInitialValueSolver <: AbstractSolver end

"""
    AbstractInitialValueSolution <: AbstractSolution

Abstract type for computed solutions of [`AbstractInitialValueProblem`](@ref)s.
"""
abstract type AbstractInitialValueSolution <: AbstractSolution end

"""
    AbstractInitialValueParameters <: AbstractParameters

Abstract type for parameter containers in [`AbstractInitialValueSolver`](@ref)s.
"""
abstract type AbstractInitialValueParameters <: AbstractParameters end

"""
    AbstractInitialValueCache <: AbstractCache

Abstract type for caching intermediate computations in [`AbstractInitialValueSolver`](@ref)s.
"""
abstract type AbstractInitialValueCache <: AbstractCache end

"""
    AbstractRightHandSide <: AbstractInitialValueParameters

Abstract type for right-hand-side (RHS) functions of [`AbstractInitialValueProblem`](@ref)s.
"""
abstract type AbstractRightHandSide <: AbstractInitialValueParameters end
