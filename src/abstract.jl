"An abstract type for all the types of [NSDE.jl](https://github.com/giancarloantonucci/NSDE.jl)."
abstract type AbstractNSDEType end

"An abstract type for all the problems of [NSDE.jl](https://github.com/giancarloantonucci/NSDE.jl)."
abstract type AbstractNSDEProblem <: AbstractNSDEType end

"An abstract type for all the solvers of [`AbstractNSDEProblem`](@ref)s."
abstract type AbstractNSDESolver <: AbstractNSDEType end

"An abstract type for the solution of an [`AbstractNSDEProblem`](@ref) obtained with an [`AbstractNSDESolver`](@ref)."
abstract type AbstractNSDESolution <: AbstractNSDEType end

"An abstract type for initial value problems."
abstract type AbstractInitialValueProblem <: AbstractNSDEProblem end

"An abstract type for solvers of [`AbstractInitialValueProblem`](@ref)s."
abstract type AbstractInitialValueSolver <: AbstractNSDESolver end

"An abstract type for the solution of an [`AbstractInitialValueProblem`](@ref) obtained with an [`AbstractInitialValueSolver`](@ref)."
abstract type AbstractInitialValueSolution <: AbstractNSDESolution end

"An abstract type for the right-hand side of an [`AbstractInitialValueProblem`](@ref)."
abstract type AbstractRightHandSideFunction <: AbstractNSDEType end
