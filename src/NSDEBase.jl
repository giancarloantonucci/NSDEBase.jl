# NSDEBase/src/NSDEBase.jl

module NSDEBase

using ArrowMacros
using LinearAlgebra
using ForwardDiff
using FiniteDifferences
using RecipesBase

include("abstract.jl")
include("utils.jl")
include("show.jl")
include("rhs/nonlinear.jl")
include("rhs/linear.jl")
include("rhs/split.jl")
include("ivp.jl")
include("odes.jl")
include("plots_recipes.jl")

"""
    solve(problem, solver; kwargs...) :: AbstractSolution

computes the solution of `problem` using `solver`, allocating a fresh cache and
solution. Solver packages extend this function; hot loops should prefer the
cache-reusing [`solve!`](@ref).
"""
function solve end

"""
    solve!(cache, solution, problem, solver) :: AbstractSolution
    solve!(solution, problem, solver) :: AbstractSolution

computes the solution of `problem` in-place, reusing `solution` and, in the
4-argument form, a pre-allocated `cache`. This is the form every hot caller
must use. Solver packages extend this function.
"""
function solve! end

"""
    initialize_cache(problem, solver) :: AbstractCache

builds a reusable cache for solving `problem` with `solver`. Solver packages
extend this function.
"""
function initialize_cache end

"""
    initialize_solution(problem, solver; kwargs...) :: AbstractSolution

builds an empty solution object for solving `problem` with `solver`. Solver
packages extend this function.
"""
function initialize_solution end

export AbstractObject

export AbstractProblem
export AbstractSolver
export AbstractSolution
export AbstractParameters
export AbstractCache

export AbstractInitialValueProblem
export AbstractInitialValueSolver
export AbstractInitialValueSolution
export AbstractInitialValueParameters
export AbstractInitialValueCache

export AbstractRightHandSide

export InitialValueProblem, IVP
export NonlinearRightHandSide, RHS
export LinearRightHandSide, LRHS
export SplitRightHandSide, SRHS

export Dahlquist
export Logistic
export SimplePendulum
export DoublePendulum
export VanDerPol
export Rössler
export Lorenz
export Lorenz96

export zero!
export copy
export solve, solve!
export initialize_cache, initialize_solution

end
