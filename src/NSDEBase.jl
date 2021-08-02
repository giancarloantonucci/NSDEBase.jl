module NSDEBase

using ArrowMacros
using ForwardDiff

abstract type NSDEProblem end
abstract type NSDESolver end
abstract type NSDESolution end

abstract type InitialValueSolver <: NSDESolver end
abstract type InitialValueSolution <: NSDESolution end

function solve end
function solve! end

include("rhs.jl")
include("ivp.jl")
include("odes.jl")

export NSDEProblem, NSDESolver, NSDESolution
export RightHandSideFunction, RHS
export InitialValueProblem, IVP, InitialValueSolver, InitialValueSolution
export Dahlquist, Logistic, Riccati, SimplePendulum, VanderPol, Rössler, Lorenz, Lorenz96
export solve, solve!

end
