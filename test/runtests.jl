using NSDEBase
using Test

@testset "IVP and RHS" begin
    f(u, t) = @. u * (1.0 - u)
    problem = IVP(f, 0.5, (0.0, 1.0))
    @test problem isa InitialValueProblem
    @test showable(MIME("text/plain"), problem)
    @test problem.rhs isa RightHandSideFunction
    @test showable(MIME("text/plain"), problem.rhs)

    f!(du, u, t) = @. du = u * (1.0 - u)
    problem = IVP(f!, [0.5], 0.0, 1.0)
    @test problem isa InitialValueProblem
    @test problem.rhs isa RightHandSideFunction
end

@testset "ODEs" begin
    @test Dahlquist() isa InitialValueProblem
    @test Logistic() isa InitialValueProblem
    @test Riccati() isa InitialValueProblem
    @test SimplePendulum() isa InitialValueProblem
    @test VanderPol() isa InitialValueProblem
    @test Rössler() isa InitialValueProblem
    @test Lorenz() isa InitialValueProblem
    @test Lorenz96() isa InitialValueProblem
end
