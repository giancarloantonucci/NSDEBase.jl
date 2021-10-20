module NSDEBase

############################################################################################
#                                         EXPORTS                                          #
############################################################################################

export AbstractNSDEProblem
export AbstractNSDESolution
export AbstractNSDESolver


export AbstractInitialValueProblem
export AbstractInitialValueSolution
export AbstractInitialValueSolver
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

############################################################################################
#                                           CORE                                           #
############################################################################################

using ArrowMacros
using LinearAlgebra
using ForwardDiff
using FiniteDifferences

include("abstract.jl")
include("utils.jl")
include("rhs.jl")
include("srhs.jl")
include("ivp.jl")
include("odes.jl")

end
