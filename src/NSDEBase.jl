module NSDEBase

export NSDEProblem, NSDESolver, NSDESolution
export RightHandSideFunction, RHS
export InitialValueProblem, IVP, InitialValueSolver, InitialValueSolution
export Dahlquist, Logistic, SimplePendulum, VanderPol, Rössler, Lorenz, Lorenz96
export solve, solve!

using ArrowMacros
using LinearAlgebra
using ForwardDiff
using FiniteDifferences

"An abstract type for the problems from [NSDE.jl](https://github.com/antonuccig/NSDE.jl)."
abstract type NSDEProblem end

"An abstract type for the solvers from [NSDE.jl](https://github.com/antonuccig/NSDE.jl)."
abstract type NSDESolver end

"An abstract type for the solutions from [NSDE.jl](https://github.com/antonuccig/NSDE.jl)."
abstract type NSDESolution end

include("utils.jl")
include("rhs.jl")
include("srhs.jl")
include("ivp.jl")
include("odes.jl")

end
