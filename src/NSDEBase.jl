module NSDEBase

using ArrowMacros
using LinearAlgebra
using ForwardDiff
using FiniteDifferences

function solve end
function solve! end

include("abstract.jl")
include("object.jl")
include("rhs.jl")
include("srhs.jl")
include("ivp.jl")
include("odes.jl")

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

end
