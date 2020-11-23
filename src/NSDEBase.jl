module NSDEBase

export BlockVector, zero!, flatview!, nestedview!
export AbstractProblem, AbstractSolver, AbstractSolution
export InitialValueProblem, IVP, RightHandSideFunction, RHS, InitialValueSolver, InitialValueSolution
export Dahlquist, Logistic, Riccati, SimplePendulum, VanderPol, Lorenz, Rössler
export solve

using ArrowMacros
using ForwardDiff: jacobian, jacobian!

abstract type NSDEProblem end
abstract type NSDESolver end
abstract type NSDESolution end

abstract type InitialValueSolver <: NSDESolver end
abstract type InitialValueSolution <: NSDESolution end

function solve end

include("arrays.jl")
include("rhs.jl")
include("ivp.jl")
include("odes.jl")

end
