abstract type AbstractShowable end

abstract type AbstractProblem <: AbstractShowable end
abstract type AbstractSolver <: AbstractShowable end
abstract type AbstractSolution <: AbstractShowable end
abstract type AbstractParameters <: AbstractShowable end
abstract type AbstractCache <: AbstractShowable end

abstract type AbstractInitialValueProblem <: AbstractProblem end
abstract type AbstractInitialValueSolver <: AbstractSolver end
abstract type AbstractInitialValueSolution <: AbstractSolution end
abstract type AbstractInitialValueParameters <: AbstractParameters end
abstract type AbstractInitialValueCache <: AbstractCache end

abstract type AbstractRightHandSide <: AbstractInitialValueParameters end
