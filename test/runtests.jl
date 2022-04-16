using NSDEBase
using Test

@testset "IVP" begin
    f(u, t) = u
    u0 = [0.0]
    tspan = (0.0, 1.0)
    problem = IVP(f, u0, tspan)
    @test problem isa InitialValueProblem
    @test problem.u0 isa AbstractVector{<:AbstractFloat}
    @test problem.tspan isa Tuple{AbstractFloat, AbstractFloat}
    # @test repr(problem) == ???
    
    f!(du, u, t) = du .= u
    u0 = t0 = tN = 0.0
    problem = IVP(f!, u0, t0, tN)
    @test problem isa InitialValueProblem
    @test problem.rhs isa NonlinearRightHandSide
    @test problem.u0 isa AbstractVector{<:AbstractFloat}
    @test problem.tspan isa Tuple{AbstractFloat, AbstractFloat}
    
    problem2 = copy(problem)
    @test problem2.rhs == problem.rhs
    @test problem2.u0 == problem.u0
    @test problem2.tspan == problem.tspan
    
    u0 = [1.0]
    t0 = tN = 1.0
    problem3 = makesub(problem, u0, t0, tN)
    @test problem3.rhs == problem.rhs
    @test problem3.u0 == u0
    @test problem3.tspan == (t0, tN)
end

@testset "RHS" begin
    f(u, t) = @. u * (1.0 - u)
    u0 = [0.0]
    tspan = (0.0, 1.0)
    problem = IVP(f, u0, tspan)
    @test problem.rhs isa NonlinearRightHandSide
    @test problem.rhs.f isa Function
    @test problem.rhs.f! isa Function
    @test problem.rhs.Df isa Function
    @test problem.rhs.Df! isa Function
    # @test repr(problem.rhs) == ???
    
    rhs2 = copy(problem.rhs)
    @test rhs2.f == problem.rhs.f
    @test rhs2.f! == problem.rhs.f!
    @test rhs2.Df == problem.rhs.Df
    @test rhs2.Df! == problem.rhs.Df!
end

@testset "ODEs" begin
    @test Dahlquist() isa InitialValueProblem
    @test Logistic() isa InitialValueProblem
    @test SimplePendulum() isa InitialValueProblem
    @test DoublePendulum() isa InitialValueProblem
    @test VanderPol() isa InitialValueProblem
    @test Rössler() isa InitialValueProblem
    @test Lorenz() isa InitialValueProblem
    @test Lorenz96() isa InitialValueProblem
end
