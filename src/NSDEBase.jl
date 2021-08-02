module NSDEBase

export NSDEProblem, NSDESolver, NSDESolution
export RightHandSideFunction, RHS
export InitialValueProblem, IVP, InitialValueSolver, InitialValueSolution
export Dahlquist, Logistic, Riccati, SimplePendulum, VanderPol, Rössler, Lorenz, Lorenz96
export solve, solve!

using ArrowMacros
using ForwardDiff

"An abstract type for problems."
abstract type NSDEProblem end

"An abstract type for solvers."
abstract type NSDESolver end

"An abstract type for solutions."
abstract type NSDESolution end

include("utils.jl")
include("rhs.jl")
include("ivp.jl")
include("odes.jl")

end
