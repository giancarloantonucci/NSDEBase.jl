var documenterSearchIndex = {"docs":
[{"location":"#NSDEBase.jl","page":"Home","title":"NSDEBase.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This is the documentation of NSDEBase.jl, a Julia package containing common types and utility functions shared by NSDE.jl and its subpackages.","category":"page"},{"location":"#API","page":"Home","title":"API","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"All exported types and functions are considered part of the public API, and thus documented in this manual.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [NSDEBase]","category":"page"},{"location":"#NSDEBase.InitialValueProblem","page":"Home","title":"NSDEBase.InitialValueProblem","text":"InitialValueProblem <: AbstractInitialValueProblem\n\nA composite type for an initial-value problem.\n\nConstructors\n\nInitialValueProblem(rhs, u0, tspan)\nInitialValueProblem(rhs, u0, t0, tN)\nIVP(args...; kwargs...)\n\nArguments\n\nrhs::AbstractRightHandSide : right-hand side function\nu0::AbstractVector : initial condition\ntspan::Tuple : end limits of time domain\n\nFunctions\n\ncopy : copy instance of problem with updated parameters\n\n\n\n\n\n","category":"type"},{"location":"#NSDEBase.LinearRightHandSide","page":"Home","title":"NSDEBase.LinearRightHandSide","text":"LinearRightHandSide <: AbstractRightHandSide\n\nA composite type for the right-hand side of an InitialValueProblem in the form f(u t) = L u + g(t).\n\nConstructors\n\nLinearRightHandSide(L, g, g!)\nLinearRightHandSide(L[, g!_or_g])\nLRHS(args...; kwargs...)\n\nArguments\n\nL::AbstractMatrix : L, the coefficient term\ng::Function : g, the forcing term, independent of u\ng!::Function : g (in-place)\n\n\n\n\n\n","category":"type"},{"location":"#NSDEBase.LinearRightHandSide-Tuple{Any, Any}","page":"Home","title":"NSDEBase.LinearRightHandSide","text":"(rhs::LinearRightHandSide)(u, t)\n(rhs::LinearRightHandSide)(du, u, t)\n\nreturns the derivative du from the solution u and time t.\n\n\n\n\n\n","category":"method"},{"location":"#NSDEBase.NonlinearRightHandSide","page":"Home","title":"NSDEBase.NonlinearRightHandSide","text":"NonlinearRightHandSide <: AbstractRightHandSide\n\nA composite type for the right-hand side of an InitialValueProblem in the generic form f(u t).\n\nConstructors\n\nNonlinearRightHandSide(f, f!, Df, Df!)\nNonlinearRightHandSide(f!_or_f, is_complex=false)\nRightHandSide(args...; kwargs...)\nRHS(args...; kwargs...)\n\nArguments\n\nf::Function : f, the right-hand side function\nf!::Function : f (in-place)\nDf::Function : mathcalDf, the Jacobian of f, with respect to u\nDf!::Function : mathcalDf (in-place)\n\n\n\n\n\n","category":"type"},{"location":"#NSDEBase.NonlinearRightHandSide-Tuple{Any, Any}","page":"Home","title":"NSDEBase.NonlinearRightHandSide","text":"(rhs::NonlinearRightHandSide)(u, t)\n(rhs::NonlinearRightHandSide)(du, u, t)\n\nreturns the derivative du from the solution u and time t.\n\n\n\n\n\n","category":"method"},{"location":"#NSDEBase.SplitRightHandSide","page":"Home","title":"NSDEBase.SplitRightHandSide","text":"SplitRightHandSide <: AbstractRightHandSide\n\nA composite type for the right-hand side of an InitialValueProblem in the form f(u t) = f_texts(u t) + f_textns(u t).\n\nConstructors\n\nSplitRightHandSide(fₛ, fₙₛ)\nSRHS(args...; kwargs...)\n\nArguments\n\nfₛ::Union{LinearRightHandSide,NonlinearRightHandSide} : f_texts, the stiff part of the right-hand side function f\nfₙₛ::NonlinearRightHandSide} : f_textns, the non-stiff part of the right-hand side function f\n\n\n\n\n\n","category":"type"},{"location":"#NSDEBase.SplitRightHandSide-Tuple{Any, Any}","page":"Home","title":"NSDEBase.SplitRightHandSide","text":"(rhs::SplitRightHandSide)(u, t)\n(rhs::SplitRightHandSide)(du, u, t)\n\nreturns the derivative du from the solution u and time t.\n\n\n\n\n\n","category":"method"},{"location":"#Base.copy-Union{Tuple{ℝ}, Tuple{ℂ}, Tuple{InitialValueProblem, Union{AbstractVector{ℂ}, ℂ}, Tuple{ℝ, ℝ}}} where {ℂ<:Number, ℝ<:Real}","page":"Home","title":"Base.copy","text":"copy(problem, u0, tspan)::InitialValueProblem\ncopy(problem, u0, t0, tN)::InitialValueProblem\n\nreturns a copy of problem with the same rhs but different u0 and tspan.\n\n\n\n\n\n","category":"method"},{"location":"#Base.show-Tuple{IO, AbstractObject}","page":"Home","title":"Base.show","text":"show(io::IO, object::AbstractObject)\n\nprints the full description of an object and its contents to a stream io.\n\n\n\n\n\n","category":"method"},{"location":"#Base.summary-Tuple{IO, AbstractObject}","page":"Home","title":"Base.summary","text":"summary(io::IO, object::AbstractObject)\n\nprints the short description of an object to a stream io.\n\n\n\n\n\n","category":"method"},{"location":"#NSDEBase.Dahlquist","page":"Home","title":"NSDEBase.Dahlquist","text":"Dahlquist(u0=0.5, tspan=(0.0, 1.0); λ=1.0)::InitialValueProblem\nDahlquist(u0, t0, tN; kwargs...)::InitialValueProblem\n\nreturns an InitialValueProblem for the Dahlquist equation.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.DoublePendulum","page":"Home","title":"NSDEBase.DoublePendulum","text":"DoublePendulum(u0=[π/4, π/4, 0.0, 0.0], tspan=(0.0, 1.0); μ=1.0, λ=1.0)::InitialValueProblem\nDoublePendulum(u0, t0, tN; kwargs...)::InitialValueProblem\n\nreturns an InitialValueProblem for the double pendulum problem.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.IVP-Tuple","page":"Home","title":"NSDEBase.IVP","text":"InitialValueProblem <: AbstractInitialValueProblem\n\nA composite type for an initial-value problem.\n\nConstructors\n\nInitialValueProblem(rhs, u0, tspan)\nInitialValueProblem(rhs, u0, t0, tN)\nIVP(args...; kwargs...)\n\nArguments\n\nrhs::AbstractRightHandSide : right-hand side function\nu0::AbstractVector : initial condition\ntspan::Tuple : end limits of time domain\n\nFunctions\n\ncopy : copy instance of problem with updated parameters\n\n\n\n\n\n","category":"method"},{"location":"#NSDEBase.LRHS-Tuple","page":"Home","title":"NSDEBase.LRHS","text":"LinearRightHandSide <: AbstractRightHandSide\n\nA composite type for the right-hand side of an InitialValueProblem in the form f(u t) = L u + g(t).\n\nConstructors\n\nLinearRightHandSide(L, g, g!)\nLinearRightHandSide(L[, g!_or_g])\nLRHS(args...; kwargs...)\n\nArguments\n\nL::AbstractMatrix : L, the coefficient term\ng::Function : g, the forcing term, independent of u\ng!::Function : g (in-place)\n\n\n\n\n\n","category":"method"},{"location":"#NSDEBase.Logistic","page":"Home","title":"NSDEBase.Logistic","text":"Logistic(u0=0.5, tspan=(0.0, 1.0); λ=1.0)::InitialValueProblem\nLogistic(u0, t0, tN; kwargs...)::InitialValueProblem\n\nreturns an InitialValueProblem for the Logistic equation.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.Lorenz","page":"Home","title":"NSDEBase.Lorenz","text":"Lorenz(u0=[2.0, 3.0, -14.0], tspan=(0.0, 1.0); σ=10.0, β=8/3, ρ=28.0)::InitialValueProblem\nLorenz(u0, t0, tN; kwargs...)::InitialValueProblem\n\nreturns an InitialValueProblem for the Lorenz equations.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.Lorenz96","page":"Home","title":"NSDEBase.Lorenz96","text":"Lorenz96(u0=[ones(39); 1.01], tspan=(0.0, 1.0); F=8.0)::InitialValueProblem\nLorenz96(u0, t0, tN; kwargs...)::InitialValueProblem\n\nreturns an InitialValueProblem for the Lorenz-96 equations.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.RHS-Tuple","page":"Home","title":"NSDEBase.RHS","text":"NonlinearRightHandSide <: AbstractRightHandSide\n\nA composite type for the right-hand side of an InitialValueProblem in the generic form f(u t).\n\nConstructors\n\nNonlinearRightHandSide(f, f!, Df, Df!)\nNonlinearRightHandSide(f!_or_f, is_complex=false)\nRightHandSide(args...; kwargs...)\nRHS(args...; kwargs...)\n\nArguments\n\nf::Function : f, the right-hand side function\nf!::Function : f (in-place)\nDf::Function : mathcalDf, the Jacobian of f, with respect to u\nDf!::Function : mathcalDf (in-place)\n\n\n\n\n\n","category":"method"},{"location":"#NSDEBase.RightHandSide-Tuple","page":"Home","title":"NSDEBase.RightHandSide","text":"NonlinearRightHandSide <: AbstractRightHandSide\n\nA composite type for the right-hand side of an InitialValueProblem in the generic form f(u t).\n\nConstructors\n\nNonlinearRightHandSide(f, f!, Df, Df!)\nNonlinearRightHandSide(f!_or_f, is_complex=false)\nRightHandSide(args...; kwargs...)\nRHS(args...; kwargs...)\n\nArguments\n\nf::Function : f, the right-hand side function\nf!::Function : f (in-place)\nDf::Function : mathcalDf, the Jacobian of f, with respect to u\nDf!::Function : mathcalDf (in-place)\n\n\n\n\n\n","category":"method"},{"location":"#NSDEBase.Rössler","page":"Home","title":"NSDEBase.Rössler","text":"Rössler(u0=[2.0, 0.0, 0.0], tspan=(0.0, 1.0); α=0.2, β=0.2, γ=5.7)::InitialValueProblem\nRössler(u0, t0, tN; kwargs...)::InitialValueProblem\n\nreturns an InitialValueProblem for the Rössler equations.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.SRHS-Tuple","page":"Home","title":"NSDEBase.SRHS","text":"SplitRightHandSide <: AbstractRightHandSide\n\nA composite type for the right-hand side of an InitialValueProblem in the form f(u t) = f_texts(u t) + f_textns(u t).\n\nConstructors\n\nSplitRightHandSide(fₛ, fₙₛ)\nSRHS(args...; kwargs...)\n\nArguments\n\nfₛ::Union{LinearRightHandSide,NonlinearRightHandSide} : f_texts, the stiff part of the right-hand side function f\nfₙₛ::NonlinearRightHandSide} : f_textns, the non-stiff part of the right-hand side function f\n\n\n\n\n\n","category":"method"},{"location":"#NSDEBase.SimplePendulum","page":"Home","title":"NSDEBase.SimplePendulum","text":"SimplePendulum(u0=[π/4, 0.0], tspan=(0.0, 2π/3 * √9.81))::InitialValueProblem\nSimplePendulum(u0, t0, tN; kwargs...)::InitialValueProblem\n\nreturns an InitialValueProblem for the simple pendulum problem.\n\n\n\n\n\n","category":"function"},{"location":"#NSDEBase.VanDerPol","page":"Home","title":"NSDEBase.VanDerPol","text":"VanDerPol(u0=[1.0, 0.0], tspan=(0.0, 1.0); μ=1.0)::InitialValueProblem\nVanDerPol(u0, t0, tN; kwargs...)::InitialValueProblem\n\nreturns an InitialValueProblem for the Van der Pol equation (in first-order form).\n\n\n\n\n\n","category":"function"}]
}
