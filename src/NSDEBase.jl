module NSDEBase

using ArrowMacros
using ForwardDiff

"Abstract type for problems solved with NSDE.jl"
abstract type NSDEProblem end
"Abstract type for solvers from NSDE.jl"
abstract type NSDESolver end
"Abstract type for solutions of problems solved with NSDE.jl"
abstract type NSDESolution end

"Abstract type for solvers of initial value problems from NSDE.jl"
abstract type InitialValueSolver <: NSDESolver end
"Abstract type for solutions of initial value problems from NSDE.jl"
abstract type InitialValueSolution <: NSDESolution end

function solve end
function solve! end

include("show.jl")
include("rhs.jl")
include("ivp.jl")
include("odes.jl")

export NSDEProblem, NSDESolver, NSDESolution
export RightHandSideFunction, RHS
export InitialValueProblem, IVP, InitialValueSolver, InitialValueSolution
export Dahlquist, Logistic, Riccati, SimplePendulum, VanderPol, Rössler, Lorenz, Lorenz96
export solve, solve!

end
