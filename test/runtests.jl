using NSDEBase
using Test

@testset "IVP and RHS" begin
    f(u, t) = @. u * (1 - u)
    u0 = 0.5
    tspan = (0, 1)
    problem = IVP(f, u0, tspan)
    @test problem isa InitialValueProblem
    @test problem.rhs isa RightHandSideFunction

    u0 = [1]
    t0 = 0
    tN = 1
    f!(du, u, t) = @. du = u * (1 - u)
    problem = IVP(f!, u0, t0, tN)
    @test problem isa InitialValueProblem
    @test problem.rhs isa RightHandSideFunction
    
    @test showable(MIME("text/plain"), problem)
    @test showable(MIME("text/plain"), problem.rhs)
end

@testset "ODEs" begin
    @test Dahlquist(0.5, 0, 1) isa InitialValueProblem
    @test Logistic(0.5, 0, 1) isa InitialValueProblem
    @test SimplePendulum([π/4, 0], 0, 1) isa InitialValueProblem
    @test DoublePendulum([π/4, π/4, 0, 0], 0, 1) isa InitialValueProblem
    @test VanderPol([1, 0], 0, 1) isa InitialValueProblem
    @test Rössler([2, 0, 0], 0, 1) isa InitialValueProblem
    @test Lorenz([2, 3, -14], 0, 1) isa InitialValueProblem
    @test Lorenz96(ones(4), 0, 1) isa InitialValueProblem
end

@testset "Functions" begin
    problem = Dahlquist()
    
    problem2 = copy(problem)
    @test problem2.rhs == problem.rhs
    @test problem2.u0 == problem.u0
    @test problem2.tspan == problem.tspan
    
    u0 = [0.0]
    t0 = tN = 0.0
    problem3 = subproblemof(problem, u0, t0, tN)
    @test problem3.rhs == problem.rhs
    @test problem3.u0 == u0
    @test problem3.tspan == (t0, tN)
end
