"An abstract type for any object from [NSDE.jl](https://github.com/antonuccig/NSDE.jl)."
abstract type AbstractNSDEObject end

"An abstract type for the problems from [NSDE.jl](https://github.com/antonuccig/NSDE.jl)."
abstract type AbstractNSDEProblem <: AbstractNSDEObject end

"An abstract type for the solvers from [NSDE.jl](https://github.com/antonuccig/NSDE.jl)."
abstract type AbstractNSDESolver <: AbstractNSDEObject end

"An abstract type for the solutions from [NSDE.jl](https://github.com/antonuccig/NSDE.jl)."
abstract type AbstractNSDESolution <: AbstractNSDEObject end

"An abstract type for initial value problems."
abstract type AbstractInitialValueProblem <: AbstractNSDEProblem end

"An abstract type for the solvers of [`AbstractInitialValueProblem`](@ref)s."
abstract type AbstractInitialValueSolver <: AbstractNSDESolver end

"An abstract type for the solution of an [`AbstractInitialValueProblem`](@ref) obtained with an [`AbstractInitialValueSolver`](@ref)."
abstract type AbstractInitialValueSolution <: AbstractNSDESolution end

"An abstract type for the right-hand side of an [`AbstractInitialValueProblem`](@ref)."
abstract type AbstractRightHandSideFunction <: AbstractNSDEObject end