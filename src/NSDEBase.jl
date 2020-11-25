module NSDEBase

export AbstractProblem, AbstractSolver, AbstractSolution
export InitialValueProblem, IVP, RightHandSideFunction, RHS, InitialValueSolver, InitialValueSolution
export Dahlquist, Logistic, Riccati, SimplePendulum, VanderPol, Lorenz, Rössler
export solve

using ForwardDiff: jacobian, jacobian!

abstract type NSDEProblem end
abstract type NSDESolver end
abstract type NSDESolution end

abstract type InitialValueSolver <: NSDESolver end
abstract type InitialValueSolution <: NSDESolution end

function solve end

include("rhs.jl")
include("ivp.jl")
include("odes.jl")

end
