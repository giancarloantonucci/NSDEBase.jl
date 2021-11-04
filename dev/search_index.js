var documenterSearchIndex = {"docs":
[{"location":"#NSDEBase.jl","page":"Home","title":"NSDEBase.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This is the documentation of NSDEBase.jl, a Julia package containing shared types and functions of NSDE.jl and its sub-packages.","category":"page"},{"location":"#Manual","page":"Home","title":"Manual","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Depth = 3","category":"page"},{"location":"#API","page":"Home","title":"API","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"All exported types and functions are considered part of the public API and thus documented in this manual.","category":"page"},{"location":"#Abstract-types","page":"Home","title":"Abstract types","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"AbstractNSDEObject\nAbstractNSDEProblem\nAbstractNSDESolver\nAbstractNSDESolution\nAbstractInitialValueProblem\nAbstractInitialValueSolver\nAbstractInitialValueSolution\nAbstractRightHandSideFunction","category":"page"},{"location":"#NSDEBase.AbstractNSDEObject","page":"Home","title":"NSDEBase.AbstractNSDEObject","text":"An abstract type for all objects from NSDE.jl.\n\n\n\n\n\n","category":"type"},{"location":"#NSDEBase.AbstractNSDEProblem","page":"Home","title":"NSDEBase.AbstractNSDEProblem","text":"An abstract type for problems from NSDE.jl.\n\n\n\n\n\n","category":"type"},{"location":"#NSDEBase.AbstractNSDESolver","page":"Home","title":"NSDEBase.AbstractNSDESolver","text":"An abstract type for solvers of AbstractNSDEProblems.\n\n\n\n\n\n","category":"type"},{"location":"#NSDEBase.AbstractNSDESolution","page":"Home","title":"NSDEBase.AbstractNSDESolution","text":"An abstract type for the solution of an AbstractNSDEProblem obtained with an AbstractNSDESolver.\n\n\n\n\n\n","category":"type"},{"location":"#NSDEBase.AbstractInitialValueProblem","page":"Home","title":"NSDEBase.AbstractInitialValueProblem","text":"An abstract type for initial value problems.\n\n\n\n\n\n","category":"type"},{"location":"#NSDEBase.AbstractInitialValueSolver","page":"Home","title":"NSDEBase.AbstractInitialValueSolver","text":"An abstract type for solvers of AbstractInitialValueProblems.\n\n\n\n\n\n","category":"type"},{"location":"#NSDEBase.AbstractInitialValueSolution","page":"Home","title":"NSDEBase.AbstractInitialValueSolution","text":"An abstract type for the solution of an AbstractInitialValueProblem obtained with an AbstractInitialValueSolver.\n\n\n\n\n\n","category":"type"},{"location":"#NSDEBase.AbstractRightHandSideFunction","page":"Home","title":"NSDEBase.AbstractRightHandSideFunction","text":"An abstract type for the right-hand side of an AbstractInitialValueProblem.\n\n\n\n\n\n","category":"type"},{"location":"#Composite-types","page":"Home","title":"Composite types","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"InitialValueProblem\nsubproblemof\nRightHandSideFunction\nSplitRightHandSideFunction","category":"page"},{"location":"#NSDEBase.InitialValueProblem","page":"Home","title":"NSDEBase.InitialValueProblem","text":"InitialValueProblem <: AbstractInitialValueProblem\n\nA composite type for an AbstractInitialValueProblem.\n\nConstructors\n\nInitialValueProblem(rhs, u0, tspan)\nInitialValueProblem(rhs, u0, t0, tN)\nIVP(args...; kwargs...)\n\nwhere\n\nrhs   :: Union{Function, AbstractRightHandSideFunction} : right-hand side derivative.\nu0    :: Union{Number, AbstractVector{<:Number}}        : initial condition.\ntspan :: Tuple{Real, Real}                              : time domain.\nt0    :: Real                                           : initial time.\ntN    :: Real                                           : final time.\n\nFunctions\n\nsubproblemof : creates a subproblem.\n\n\n\n\n\n","category":"type"},{"location":"#NSDEBase.subproblemof","page":"Home","title":"NSDEBase.subproblemof","text":"subproblemof(problem::InitialValueProblem, u0, tspan) :: InitialValueProblem\nsubproblemof(problem, u0, t0, tN) :: InitialValueProblem\n\nreturns a copy of problem with different initial condition and time domain.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.RightHandSideFunction","page":"Home","title":"NSDEBase.RightHandSideFunction","text":"RightHandSideFunction <: AbstractRightHandSideFunction\n\nA composite type for the right-hand side of an InitialValueProblem.\n\nConstructors\n\nRightHandSideFunction(f, f!, Df, Df!)\nRightHandSideFunction(f!_or_f)\nRHS(args...; kwargs...)\n\nArguments\n\nf       :: Function : right-hand side derivative.\nf!      :: Function : right-hand side derivative (in-place).\nDf      :: Function : jacobian of right-hand side derivative.\nDf!     :: Function : jacobian of right-hand side derivative (in-place).\nf!_or_f :: Function : function (in-place or not) from which all other fields will be constructed.\n\n\n\n\n\n","category":"type"},{"location":"#NSDEBase.SplitRightHandSideFunction","page":"Home","title":"NSDEBase.SplitRightHandSideFunction","text":"SplitRightHandSideFunction <: AbstractRightHandSideFunction\n\nA composite type for the split right-hand side of an InitialValueProblem.\n\nConstructors\n\nSplitRightHandSideFunction(L, rhs)\nSRHS(args...; kwargs...)\n\nArguments\n\nL   :: Union{Number, AbstractMatrix{<:Number}} : linear part of right-hand side derivative.\nrhs :: RightHandSideFunction                   : nonlinear part of right-hand side derivative.\n\n\n\n\n\n","category":"type"},{"location":"#ODEs","page":"Home","title":"ODEs","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Dahlquist\nLogistic\nSimplePendulum\nVanderPol\nRössler\nLorenz\nLorenz96","category":"page"},{"location":"#NSDEBase.Dahlquist","page":"Home","title":"NSDEBase.Dahlquist","text":"Dahlquist(u0=0.5, tspan=(0.0, 1.0); λ=1.0) :: InitialValueProblem\nDahlquist(u0, t0, tN; kwargs...) :: InitialValueProblem\n\nreturns an InitialValueProblem for the Dahlquist equation.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.Logistic","page":"Home","title":"NSDEBase.Logistic","text":"Logistic(u0=0.5, tspan=(0.0, 1.0); λ=1.0) :: InitialValueProblem\nLogistic(u0, t0, tN; kwargs...) :: InitialValueProblem\n\nreturns an InitialValueProblem for the Logistic equation.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.SimplePendulum","page":"Home","title":"NSDEBase.SimplePendulum","text":"SimplePendulum(u0=[0.0, π/2], tspan=(0.0, 2π); g=1.0, L=1.0) :: InitialValueProblem\nSimplePendulum(u0, t0, tN; kwargs...) :: InitialValueProblem\n\nreturns an InitialValueProblem for the simple pendulum problem.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.VanderPol","page":"Home","title":"NSDEBase.VanderPol","text":"VanderPol(u0=[1.0, 0.0], tspan=(0.0, 1.0); μ=1.0) :: InitialValueProblem\nVanderPol(u0, t0, tN; kwargs...) :: InitialValueProblem\n\nreturns an InitialValueProblem for the Van der Pol equation (in first-order form).\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.Rössler","page":"Home","title":"NSDEBase.Rössler","text":"Rössler(u0=[2.0, 0.0, 0.0], tspan=(0.0, 1.0); a=0.2, b=0.2, c=5.7) :: InitialValueProblem\nRössler(u0, t0, tN; kwargs...) :: InitialValueProblem\n\nreturns an InitialValueProblem for the Rössler equations.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.Lorenz","page":"Home","title":"NSDEBase.Lorenz","text":"Lorenz(u0=[2.0, 3.0, -14.0], tspan=(0.0, 1.0); σ=10.0, β=8/3, ρ=28.0) :: InitialValueProblem\nLorenz(u0, t0, tN; kwargs...) :: InitialValueProblem\n\nreturns an InitialValueProblem for the Lorenz equations.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.Lorenz96","page":"Home","title":"NSDEBase.Lorenz96","text":"Lorenz96(u0=[1.01; ones(40)], tspan=(0.0, 1.0); F=8) :: InitialValueProblem\nLorenz96(u0, t0, tN; kwargs...) :: InitialValueProblem\n\nreturns an InitialValueProblem for the Lorenz-96 equations.\n\n\n\n\n\n","category":"function"},{"location":"#Utilities","page":"Home","title":"Utilities","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"show\nsummary","category":"page"},{"location":"#Base.show","page":"Home","title":"Base.show","text":"show(io::IO, object::AbstractNSDEObject)\n\nprints the full description of an object and its contents to a stream io.\n\n\n\n\n\n","category":"function"},{"location":"#Base.summary","page":"Home","title":"Base.summary","text":"summary(io::IO, problem::AbstractNSDEObject)\n\nprints the short description of an object to a stream io.\n\n\n\n\n\n","category":"function"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}
