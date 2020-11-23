using NSDEBase
using Test

@testset "RHS and IVP" begin
    f(u, t) = @. u * (1.0 - u)
    problem = IVP(f, 0.5, 0.0, 1.0)
    @test problem isa InitialValueProblem
end

@testset "ODEs" begin
    problem = Dahlquist()
    @test problem isa InitialValueProblem
    problem = Logistic()
    @test problem isa InitialValueProblem
    problem = Riccati()
    @test problem isa InitialValueProblem
    problem = SimplePendulum()
    @test problem isa InitialValueProblem
    problem = VanderPol()
    @test problem isa InitialValueProblem
    problem = Lorenz()
    @test problem isa InitialValueProblem
    problem = Rössler()
    @test problem isa InitialValueProblem
end
