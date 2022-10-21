var documenterSearchIndex = {"docs":
[{"location":"#NSDEBase.jl","page":"Home","title":"NSDEBase.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This is the documentation of NSDEBase.jl, a Julia package containing common types and utility functions shared by NSDE.jl and its subpackages.","category":"page"},{"location":"#Manual","page":"Home","title":"Manual","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Depth = 3","category":"page"},{"location":"#API","page":"Home","title":"API","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"All exported types and functions are considered part of the public API, and thus documented in this manual.","category":"page"},{"location":"#Abstract-types","page":"Home","title":"Abstract types","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"AbstractNSDEType\nAbstractNSDEProblem\nAbstractNSDESolver\nAbstractNSDESolution\nAbstractInitialValueProblem\nAbstractInitialValueSolver\nAbstractInitialValueSolution\nAbstractRightHandSide","category":"page"},{"location":"#Composite-types","page":"Home","title":"Composite types","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"InitialValueProblem\nsubproblemof\nRightHandSide\nSplitRightHandSide","category":"page"},{"location":"#NSDEBase.InitialValueProblem","page":"Home","title":"NSDEBase.InitialValueProblem","text":"InitialValueProblem <: AbstractInitialValueProblem\n\nA composite type for an initial-value problem.\n\nConstructors\n\nInitialValueProblem(rhs, u0, tspan)\nInitialValueProblem(rhs, u0, t0, tN)\nIVP(args...; kwargs...)\n\nArguments\n\nrhs :: Union{AbstractRightHandSide, Function, AbstractMatrix{ℂ}, ℂ} where ℂ<:Number : right-hand side function.\nu0 :: Union{AbstractVector{ℂ}, ℂ} where ℂ<:Number : initial condition.\ntspan :: Tuple{ℝ, ℝ} where ℝ<:Real : end limits of time domain.\n\nFunctions\n\nsubproblemof : creates a subproblem.\n\n\n\n\n\n","category":"type"},{"location":"#NSDEBase.subproblemof","page":"Home","title":"NSDEBase.subproblemof","text":"subproblemof(problem, u0, tspan)  :: InitialValueProblem\nsubproblemof(problem, u0, t0, tN) :: InitialValueProblem\n\nreturns a subproblem of problem with the same rhs but different u0 and tspan.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.SplitRightHandSide","page":"Home","title":"NSDEBase.SplitRightHandSide","text":"SplitRightHandSide <: AbstractRightHandSide\n\nA composite type for the right-hand side of an InitialValueProblem in the form f_s(u t) + f_ns(u t).\n\nConstructors\n\nSplitRightHandSide(stiff, nonstiff)\nSRHS(args...; kwargs...)\n\nArguments\n\nstiff :: Union{NonlinearRightHandSide, LinearRightHandSide, Function, AbstractMatrix{ℂ}, ℂ} where ℂ<:Number : f_s.\nnonstiff :: Union{NonlinearRightHandSide, Function} : f_ns.\n\n\n\n\n\n","category":"type"},{"location":"#ODEs","page":"Home","title":"ODEs","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Dahlquist\nLogistic\nSimplePendulum\nVanderPol\nRössler\nLorenz\nLorenz96","category":"page"},{"location":"#NSDEBase.Dahlquist","page":"Home","title":"NSDEBase.Dahlquist","text":"Dahlquist(u0=0.5, tspan=(0.0, 1.0); λ=1.0) :: InitialValueProblem\nDahlquist(u0, t0, tN; kwargs...) :: InitialValueProblem\n\nreturns an InitialValueProblem for the Dahlquist equation.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.Logistic","page":"Home","title":"NSDEBase.Logistic","text":"Logistic(u0=0.5, tspan=(0.0, 1.0); λ=1.0) :: InitialValueProblem\nLogistic(u0, t0, tN; kwargs...) :: InitialValueProblem\n\nreturns an InitialValueProblem for the Logistic equation.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.SimplePendulum","page":"Home","title":"NSDEBase.SimplePendulum","text":"SimplePendulum(u0=[π/4, 0.0], tspan=(0.0, 2π/3 * √9.81)) :: InitialValueProblem\nSimplePendulum(u0, t0, tN; kwargs...) :: InitialValueProblem\n\nreturns an InitialValueProblem for the simple pendulum problem.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.Rössler","page":"Home","title":"NSDEBase.Rössler","text":"Rössler(u0=[2.0, 0.0, 0.0], tspan=(0.0, 1.0); α=0.2, β=0.2, γ=5.7) :: InitialValueProblem\nRössler(u0, t0, tN; kwargs...) :: InitialValueProblem\n\nreturns an InitialValueProblem for the Rössler equations.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.Lorenz","page":"Home","title":"NSDEBase.Lorenz","text":"Lorenz(u0=[2.0, 3.0, -14.0], tspan=(0.0, 1.0); σ=10.0, β=8/3, ρ=28.0) :: InitialValueProblem\nLorenz(u0, t0, tN; kwargs...) :: InitialValueProblem\n\nreturns an InitialValueProblem for the Lorenz equations.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.Lorenz96","page":"Home","title":"NSDEBase.Lorenz96","text":"Lorenz96(u0=[ones(39); 1.01], tspan=(0.0, 1.0); F=8.0) :: InitialValueProblem\nLorenz96(u0, t0, tN; kwargs...) :: InitialValueProblem\n\nreturns an InitialValueProblem for the Lorenz-96 equations.\n\n\n\n\n\n","category":"function"},{"location":"#Utilities","page":"Home","title":"Utilities","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"show\nsummary","category":"page"},{"location":"#Base.show","page":"Home","title":"Base.show","text":"show(io::IO, object::AbstractObject)\n\nprints the full description of an object and its contents to a stream io.\n\n\n\n\n\n","category":"function"},{"location":"#Base.summary","page":"Home","title":"Base.summary","text":"summary(io::IO, object::AbstractObject)\n\nprints the short description of an object to a stream io.\n\n\n\n\n\n","category":"function"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}
