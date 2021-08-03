var documenterSearchIndex = {"docs":
[{"location":"#NSDEBase.jl","page":"Home","title":"NSDEBase.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This is the documentation of NSDEBase.jl, a Julia package containing shared types and functions of NSDE.jl and its sub-packages.","category":"page"},{"location":"#Manual","page":"Home","title":"Manual","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Depth = 3","category":"page"},{"location":"#API","page":"Home","title":"API","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"All exported types and functions are considered part of the public API and thus documented in this manual.","category":"page"},{"location":"#Types","page":"Home","title":"Types","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"NSDEProblem\nNSDESolver\nNSDESolution\nRightHandSideFunction\nInitialValueProblem\nInitialValueSolver\nInitialValueSolution","category":"page"},{"location":"#NSDEBase.NSDEProblem","page":"Home","title":"NSDEBase.NSDEProblem","text":"An abstract type for the problems from NSDE.jl.\n\n\n\n\n\n","category":"type"},{"location":"#NSDEBase.NSDESolver","page":"Home","title":"NSDEBase.NSDESolver","text":"An abstract type for the solvers from NSDE.jl.\n\n\n\n\n\n","category":"type"},{"location":"#NSDEBase.NSDESolution","page":"Home","title":"NSDEBase.NSDESolution","text":"An abstract type for the solutions from NSDE.jl.\n\n\n\n\n\n","category":"type"},{"location":"#NSDEBase.RightHandSideFunction","page":"Home","title":"NSDEBase.RightHandSideFunction","text":"RightHandSideFunction\n\nA composite type for the right-hand side of an InitialValueProblem.\n\nConstructors\n\nRightHandSideFunction(f, f!, Df, Df!)\nRightHandSideFunction(f!_or_f)\nIVP(args...; kwargs...)\n\nArguments\n\nf :: Function : right-hand side derivative.\nf! :: Function : right-hand side derivative (in-place).\nDf :: Function : jacobian of right-hand side derivative.\nDf! :: Function : jacobian of right-hand side derivative (in-place).\nf!_or_f :: Function : function (in-place or not) from which all other fields are automatically constructed.\n\nFunctions\n\nshow : shows name and contents.\nsummary : shows name.\n\n\n\n\n\n","category":"type"},{"location":"#NSDEBase.InitialValueProblem","page":"Home","title":"NSDEBase.InitialValueProblem","text":"InitialValueProblem <: NSDEProblem\n\nA composite type for initial value problems.\n\nConstructors\n\nInitialValueProblem(rhs, u0, tspan)\nInitialValueProblem(rhs, u0, t0, tN)\nIVP(args...; kwargs...)\n\nArguments\n\nrhs :: Union{Function, RightHandSideFunction} : right-hand side derivative.\nu0 :: Union{Number, AbstractVector{Number}} : initial condition.\ntspan :: Tuple{Real, Real} : time domain.\nt0 :: Real : initial time.\ntN :: Real : final time.\n\nFunctions\n\nshow : shows name and contents.\nsummary : shows name.\n\n\n\n\n\n","category":"type"},{"location":"#NSDEBase.InitialValueSolver","page":"Home","title":"NSDEBase.InitialValueSolver","text":"An abstract type for the solvers of InitialValueProblems.\n\n\n\n\n\n","category":"type"},{"location":"#NSDEBase.InitialValueSolution","page":"Home","title":"NSDEBase.InitialValueSolution","text":"An abstract type for the solution of an InitialValueProblem.\n\n\n\n\n\n","category":"type"},{"location":"#ODEs","page":"Home","title":"ODEs","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Dahlquist\nLogistic\nRiccati\nSimplePendulum\nVanderPol\nRössler\nLorenz\nLorenz96","category":"page"},{"location":"#NSDEBase.Dahlquist","page":"Home","title":"NSDEBase.Dahlquist","text":"Dahlquist(u0 = 1.0, tspan = (0.0, 1.0); λ = 1.0) :: InitialValueProblem\n\nreturns an InitialValueProblem for the Dahlquist equation.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.Logistic","page":"Home","title":"NSDEBase.Logistic","text":"Logistic(u0 = 0.5, tspan = (0.0, 5.0); λ = 1.0) :: InitialValueProblem\n\nreturns an InitialValueProblem for the Logistic equation.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.Riccati","page":"Home","title":"NSDEBase.Riccati","text":"Riccati(u0 = 0.0, tspan = (0.0, 5.0)) :: InitialValueProblem\n\nreturns an InitialValueProblem for the Riccati equation.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.SimplePendulum","page":"Home","title":"NSDEBase.SimplePendulum","text":"SimplePendulum(u0 = [0.0, π/2], tspan = (0.0, 2π); g = 1.0, L = 1.0) :: InitialValueProblem\n\nreturns an InitialValueProblem for the simple pendulum problem.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.VanderPol","page":"Home","title":"NSDEBase.VanderPol","text":"VanderPol(u0 = [1.0; 0.0], tspan = (0.0, 5.0); μ = 1.0) :: InitialValueProblem\n\nreturns an InitialValueProblem for the Van der Pol equation (in first-order form).\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.Rössler","page":"Home","title":"NSDEBase.Rössler","text":"Rössler(u0 = [2.0, 0.0, 0.0], tspan = (0.0, 10.0); a = 0.2, b = 0.2, c = 5.7) :: InitialValueProblem\n\nreturns an InitialValueProblem for the Rössler equations.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.Lorenz","page":"Home","title":"NSDEBase.Lorenz","text":"Lorenz(u0 = [2.0, 3.0, -14.0], tspan = (0.0, 10.0); σ = 10.0, β = 8/3, ρ = 28.0) :: InitialValueProblem\n\nreturns an InitialValueProblem for the Lorenz equations.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.Lorenz96","page":"Home","title":"NSDEBase.Lorenz96","text":"Lorenz96(u0, tspan; F = 8) :: InitialValueProblem\n\nreturns an InitialValueProblem for the Rössler equations.\n\n\n\n\n\n","category":"function"},{"location":"#Utilities","page":"Home","title":"Utilities","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"show\nsummary","category":"page"},{"location":"#Base.show","page":"Home","title":"Base.show","text":"show(io::IO, rhs::RightHandSideFunction)\n\nprints a full description of rhs and its contents to a stream io.\n\n\n\n\n\nshow(io::IO, problem::InitialValueProblem)\n\nprints a full description of problem and its contents to a stream io.\n\n\n\n\n\n","category":"function"},{"location":"#Base.summary","page":"Home","title":"Base.summary","text":"summary(io::IO, rhs::RightHandSideFunction)\n\nprints a brief description of rhs to a stream io.\n\n\n\n\n\nsummary(io::IO, problem::InitialValueProblem)\n\nprints a brief description of problem to a stream io.\n\n\n\n\n\n","category":"function"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}
