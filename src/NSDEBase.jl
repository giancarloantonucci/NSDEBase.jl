module NSDEBase

using ArrowMacros
using LinearAlgebra
using ForwardDiff
using FiniteDifferences
using RecipesBase

function solve end
function solve! end

include("abstract.jl")
include("show.jl")
include("rhs.jl")
include("srhs.jl")
include("ivp.jl")
include("odes.jl")
include("plotrecipes.jl")

export AbstractNSDEObject
export AbstractNSDEProblem
export AbstractNSDESolver
export AbstractNSDESolution
export AbstractInitialValueProblem
export AbstractInitialValueSolver
export AbstractInitialValueSolution
export AbstractRightHandSideFunction

export InitialValueProblem, IVP
export RightHandSideFunction, RHS
export SplitRightHandSideFunction, SRHS

export Dahlquist
export Logistic
export SimplePendulum
export VanderPol
export Rössler
export Lorenz
export Lorenz96

export solve, solve!
export subproblemof

end
