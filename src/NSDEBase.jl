module NSDEBase

export NSDEProblem, NSDESolver, NSDESolution, NSDEParameters
export RightHandSideFunction, RHS
export InitialValueProblem, IVP, InitialValueSolver, InitialValueSolution, InitialValueParameters
export Dahlquist, Logistic, Riccati, SimplePendulum, VanderPol, Lorenz, Rössler
export solve

using ArrowMacros
using ForwardDiff
using Suppressor: @suppress_err

abstract type NSDEProblem end
abstract type NSDESolver end
abstract type NSDESolution end
abstract type NSDEParameters end

abstract type InitialValueSolver <: NSDESolver end
abstract type InitialValueSolution <: NSDESolution end
abstract type InitialValueParameters <: NSDEParameters end

function solve end

include("rhs.jl")
include("ivp.jl")
include("odes.jl")

end
