abstract type AbstractObject end

abstract type AbstractProblem <: AbstractObject end
abstract type AbstractSolver <: AbstractObject end
abstract type AbstractSolution <: AbstractObject end
abstract type AbstractParameters <: AbstractObject end
abstract type AbstractCache <: AbstractObject end

abstract type AbstractInitialValueProblem <: AbstractProblem end
abstract type AbstractInitialValueSolver <: AbstractSolver end
abstract type AbstractInitialValueSolution <: AbstractSolution end
abstract type AbstractInitialValueParameters <: AbstractParameters end
abstract type AbstractInitialValueCache <: AbstractCache end

abstract type AbstractRightHandSide <: AbstractInitialValueParameters end
