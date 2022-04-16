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
include("plot.jl")

function solve end
function solve! end

export AbstractShowable

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
export NonlinearRightHandSide, NRHS, RHS
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
export subproblemof
export solve, solve!

end
