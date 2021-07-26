using NSDEBase
using Test

@testset "IVP and RHS" begin
    f(u, t) = @. u * (1.0 - u)
    problem = IVP(f, 0.5, (0.0, 1.0))
    @test problem isa InitialValueProblem
    @test problem.rhs isa RightHandSideFunction

    f!(du, u, t) = @. du = u * (1.0 - u)
    problem = IVP(f!, [0.5], 0.0, 1.0)
    @test problem isa InitialValueProblem
    @test problem.rhs isa RightHandSideFunction

    problem₂ = copy(problem)
    @test problem₂.rhs == problem.rhs
    @test problem₂.u0 == problem.u0
    @test problem₂.tspan == problem.tspan
end

@testset "ODEs" begin
    @test Dahlquist() isa InitialValueProblem
    @test Logistic() isa InitialValueProblem
    @test Riccati() isa InitialValueProblem
    @test SimplePendulum() isa InitialValueProblem
    @test VanderPol() isa InitialValueProblem
    @test Lorenz() isa InitialValueProblem
    @test Rössler() isa InitialValueProblem
end
