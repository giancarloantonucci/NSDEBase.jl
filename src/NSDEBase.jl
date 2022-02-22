module NSDEBase

using ArrowMacros
using LinearAlgebra
using ForwardDiff
using FiniteDifferences
using RecipesBase

include("abstract.jl")
include("utils.jl")
include("show.jl")
include("rhs.jl")
include("srhs.jl")
include("ivp.jl")
include("odes.jl")
include("plotrecipes.jl")

export AbstractNSDEType
export AbstractNSDEProblem
export AbstractNSDESolver
export AbstractNSDESolution
export AbstractInitialValueProblem
export AbstractInitialValueSolver
export AbstractInitialValueSolution
export AbstractRightHandSideFunction

export InitialValueProblem, IVP, makesub
export RightHandSideFunction, RHS
export SplitRightHandSideFunction, SRHS

export Dahlquist
export Logistic
export SimplePendulum
export DoublePendulum
export VanderPol
export Rössler
export Lorenz
export Lorenz96

export solve, solve!

end
