using NSDEBase
using LinearAlgebra
using Test
using Aqua

# Measure allocations of an RHS call from function scope (a bare `@allocated`
# at testset level picks up global-scope dispatch noise). Call twice, keep the
# second measurement so compilation is not counted.
function measure_rhs_allocs(rhs, du, u, t)
    rhs(du, u, t)
    return @allocated rhs(du, u, t)
end
function measure_rhs_allocs(rhs, du, v, u, t)
    rhs(du, v, u, t)
    return @allocated rhs(du, v, u, t)
end
function measure_copy_allocs(problem, u0, tspan)
    copy(problem, u0, tspan)
    return @allocated copy(problem, u0, tspan)
end

@testset "Aqua" begin
    # `RecipesBase.recipetype(::Val{:phaseplot}, …)` and the :convergence twin
    # are technically piracy (no owned type in the signature) but are the
    # established idiom for user plots; each Val symbol is namespaced by our
    # own plot name, so a clash is not a real-world risk. Whitelist that one
    # function and keep every other piracy check strict.
    Aqua.test_all(NSDEBase; piracies=(; treat_as_own=[NSDEBase.RecipesBase.recipetype]))
end

@testset "NonlinearRightHandSide" begin
    @testset "out-of-place entry" begin
        f(u, t) = @. u * (1.0 - u)
        rhs = RHS(f)
        u = [0.3, -0.2]
        t = 0.7
        du_ref = [0.3 * 0.7, -0.2 * 1.2]           # u(1-u) by hand
        @test rhs(u, t) ≈ du_ref
        du = similar(u)
        @test rhs(du, u, t) ≈ du_ref               # f! derived from f
        @test du ≈ du_ref
        J_ref = Diagonal(1.0 .- 2.0 .* u)          # d/du [u(1-u)] = 1 - 2u
        @test rhs.Df(u, t) ≈ J_ref
        J = zeros(2, 2)
        rhs.Df!(J, similar(u), u, t)
        @test J ≈ J_ref
    end

    @testset "in-place entry" begin
        function f!(du, u, t)
            du[1] = u[1] * u[2]
            du[2] = sin(u[1]) + t
            return du
        end
        rhs = RHS(f!)
        u = [0.4, 1.5]
        t = 0.25
        du_ref = [0.6, sin(0.4) + 0.25]
        du = similar(u)
        @test rhs(du, u, t) ≈ du_ref
        @test rhs(u, t) ≈ du_ref                   # f derived from f!
        J_ref = [u[2] u[1]; cos(u[1]) 0.0]         # Jacobian by hand
        @test rhs.Df(u, t) ≈ J_ref
        J = zeros(2, 2)
        rhs.Df!(J, similar(u), u, t)
        @test J ≈ J_ref
        # 4-arg scratch form is a pass-through for a nonlinear RHS
        v = similar(u)
        @test rhs(du, v, u, t) ≈ du_ref
        @test measure_rhs_allocs(rhs, du, v, u, t) == 0
    end

    @testset "in-place entry that returns nothing" begin
        # A mutating function owes nothing to its caller but the filled buffer.
        # The derived out-of-place `f` used to hand back `f!`'s return value,
        # so this `rhs(u, t)` returned `nothing` and the AD Jacobian failed.
        function g!(du, u, t)
            du[1] = u[1] * u[2]
            du[2] = sin(u[1]) + t
            return nothing
        end
        rhs = RHS(g!)
        u = [0.4, 1.5]
        t = 0.25
        du_ref = [0.6, sin(0.4) + 0.25]
        @test rhs(u, t) isa AbstractVector
        @test rhs(u, t) ≈ du_ref
        @test rhs.Df(u, t) ≈ [u[2] u[1]; cos(u[1]) 0.0]
    end

    @testset "complex (Wirtinger) entry" begin
        λ = 2.0 - 1.0im
        c = 0.5 + 0.3im
        f(u, t) = @. λ * u + c * conj(u)
        rhs = RHS(f; iscomplex=true)
        u = [1.0 + 2.0im, -0.7 + 0.1im]
        t = 0.0
        @test rhs(u, t) ≈ λ .* u .+ c .* conj.(u)
        # ∂f/∂u of λu + c·conj(u) is exactly λ·I in the Wirtinger sense
        @test rhs.Df(u, t) ≈ λ * Matrix(I, 2, 2) atol = 1e-7
        # Holomorphic case: f(u) = u², ∂f/∂u = 2u on the diagonal
        rhs2 = RHS((u, t) -> u .^ 2; iscomplex=true)
        @test rhs2.Df(u, t) ≈ Diagonal(2.0 .* u) atol = 1e-7
    end

    @testset "input validation" begin
        @test_throws ArgumentError RHS(() -> 0.0)
    end
end

@testset "LinearRightHandSide" begin
    L = [0.0 1.0; -4.0 -0.4]
    u = [0.5, -1.0]
    t = 0.3

    @testset "L only" begin
        rhs = LRHS(L)
        @test rhs(u, t) ≈ L * u
        du = similar(u)
        @test rhs(du, u, t) ≈ L * u
        v = similar(u)
        @test rhs(du, v, u, t) ≈ L * u
        @test measure_rhs_allocs(rhs, du, u, t) == 0     # no g ⇒ 3-arg already clean
        @test measure_rhs_allocs(rhs, du, v, u, t) == 0
    end

    @testset "L and forcing g" begin
        g(t) = [sin(t), cos(t)]
        g!(dg, t) = (dg[1] = sin(t); dg[2] = cos(t); dg)   # allocation-free by hand
        du_ref = L * u + g(t)
        for rhs in (LRHS(L, g), LRHS(L, g!))               # g and g! entry forms
            @test rhs(u, t) ≈ du_ref
            du = similar(u)
            @test rhs(du, u, t) ≈ du_ref
            v = similar(u)
            @test rhs(du, v, u, t) ≈ du_ref
        end
        # The 4-arg form must not allocate — provable only when the user's g!
        # is itself allocation-free (a g! derived from g must call g).
        rhs = LRHS(L, g!)
        du, v = similar(u), similar(u)
        @test measure_rhs_allocs(rhs, du, v, u, t) == 0
    end

    @testset "scalar L" begin
        rhs = LRHS(2.0)
        @test rhs.L == hcat(2.0)
        @test rhs([3.0], t) ≈ [6.0]
    end

    @testset "input validation" begin
        @test_throws ArgumentError LRHS(L, () -> 0.0)
    end
end

@testset "SplitRightHandSide" begin
    L = [-2.0 0.0; 0.0 -3.0]
    g(t) = [t, 2t]
    fₙₛ(u, t) = @. sin(u)
    u = [0.9, -0.4]
    t = 0.6
    du_ref = L * u + g(t) + sin.(u)

    @testset "constructor combinations" begin
        combos = (
            SRHS(LRHS(L, g), RHS(fₙₛ)),
            SRHS(LRHS(L, g), fₙₛ),
            SRHS(RHS((u, t) -> L * u + g(t)), fₙₛ),        # nonlinear stiff part
            SRHS((u, t) -> L * u + g(t), RHS(fₙₛ)),        # Function stiff part
        )
        for rhs in combos
            @test rhs isa SplitRightHandSide
            @test rhs(u, t) ≈ du_ref
            du = similar(u)
            @test rhs(du, u, t) ≈ du_ref
            v = similar(u)
            @test rhs(du, v, u, t) ≈ du_ref
        end
        # Scalar and matrix stiff parts build an LRHS underneath
        @test SRHS(-2.0, fₙₛ).fₛ isa LinearRightHandSide
        @test SRHS(L, RHS(fₙₛ)).fₛ isa LinearRightHandSide
    end

    @testset "allocation-free 4-arg path" begin
        # The case that used to allocate twice per call: LRHS-with-g stiff part
        rhs = SRHS(LRHS(L, (dg, t) -> (dg[1] = t; dg[2] = 2t; dg)), RHS((du, u, t) -> du .= sin.(u)))
        du, v = similar(u), similar(u)
        @test rhs(du, v, u, t) ≈ du_ref
        @test measure_rhs_allocs(rhs, du, v, u, t) == 0
    end
end

@testset "InitialValueProblem" begin
    @testset "construction forms" begin
        f!(du, u, t) = du .= u
        problem = IVP(f!, 0.5, 0.0, 1.0)              # scalar u0, (t0, tN)
        @test problem isa InitialValueProblem
        @test problem.rhs isa NonlinearRightHandSide
        @test problem.u0 == [0.5]
        @test problem.tspan == (0.0, 1.0)
        @test IVP(2.0, 0.5, (0.0, 1.0)).rhs isa LinearRightHandSide
    end

    @testset "SRHS construction path" begin
        L = [-1.0 0.0; 0.0 -2.0]
        fₙₛ(u, t) = @. cos(u)
        u0 = [0.1, 0.2]
        problem = IVP(L, fₙₛ, u0, (0.0, 1.0))
        @test problem.rhs isa SplitRightHandSide
        @test problem.rhs(u0, 0.0) ≈ L * u0 + cos.(u0)
        @test IVP(L, fₙₛ, u0, 0.0, 1.0).tspan == (0.0, 1.0)
        @test IVP(-1.0, fₙₛ, 0.1, (0.0, 1.0)).u0 == [0.1]      # scalar forms
        srhs = SRHS(LRHS(L), fₙₛ)
        @test IVP(srhs, u0, (0.0, 1.0)).rhs === srhs           # direct pass-through
    end

    @testset "copy semantics" begin
        problem = Logistic([0.3], (0.0, 2.0))
        u0 = [0.9]
        sub = copy(problem, u0, 0.5, 1.5)
        @test sub.rhs === problem.rhs                 # rhs shared by identity, not rebuilt
        @test sub.u0 === u0
        @test sub.tspan == (0.5, 1.5)
        @test copy(problem, (0.5, 1.5)).u0 === problem.u0
        @test (@inferred copy(problem, u0, (0.5, 1.5))) isa typeof(problem)
        # Hot seam: building a chunk problem must stay cheap (one small struct)
        @test measure_copy_allocs(problem, u0, (0.5, 1.5)) ≤ 128
    end
end

@testset "ODEs" begin
    s = sin(π / 4)

    @testset "hand-checked derivatives at u0, t0" begin
        checks = (
            (Dahlquist(),      [0.5]),                                  # λ·u0
            (Logistic(),       [0.25]),                                 # u0(1−u0)
            (SimplePendulum(), [0.0, -s]),                              # [ω, −sin θ]
            (DoublePendulum(), [0.0, 0.0, -s, 0.0]),                    # M\v with θ₁=θ₂
            (VanDerPol(),      [0.0, -1.0]),                            # [0, −u₁]
            (Rössler(),        [0.0, 2.0, 0.2]),                        # by substitution
            (Lorenz(),         [10.0, 81.0, 6.0 + 112.0 / 3.0]),        # by substitution
        )
        for (problem, du_ref) in checks
            du = similar(problem.u0)
            @test problem.rhs(du, problem.u0, problem.tspan[1]) ≈ du_ref
        end
    end

    @testset "Lorenz96" begin
        problem = Lorenz96()                          # u0 = [1.01; ones(39)], F = 8
        du = similar(problem.u0)
        problem.rhs(du, problem.u0, 0.0)
        @test du[1] ≈ 6.99                            # (u₂−u₃₉)u₄₀ − u₁ + F
        @test du[2] ≈ 7.0                             # (u₃−u₄₀)u₁ − u₂ + F
        @test du[3] ≈ 6.99                            # (u₄−u₁)u₂ − u₃ + F
        @test du[40] ≈ 7.01                           # (u₁−u₃₈)u₃₉ − u₄₀ + F
        @test_throws ArgumentError Lorenz96(ones(3))
    end

    @testset "three-positional forms have defaults" begin
        @test Dahlquist(0.5, 0.0, 1.0) isa InitialValueProblem
        @test Logistic(0.5, 0.0, 1.0) isa InitialValueProblem
        @test DoublePendulum([π/4, π/4, 0.0, 0.0], 0.0, 1.0) isa InitialValueProblem
        @test VanDerPol([1.0, 0.0], 0.0, 1.0) isa InitialValueProblem
        @test Rössler([2.0, 0.0, 0.0], 0.0, 1.0) isa InitialValueProblem
        @test Lorenz([2.0, 3.0, -14.0], 0.0, 1.0) isa InitialValueProblem
        @test Lorenz96([1.01; ones(39)], 0.0, 1.0) isa InitialValueProblem
    end
end

@testset "utilities and show" begin
    v = [1.0, 2.0]
    @test zero!(v) == [0.0, 0.0]
    w = [[1.0], [2.0, 3.0]]
    @test zero!(w) == [[0.0], [0.0, 0.0]]

    str = sprint(show, Dahlquist())
    @test occursin("InitialValueProblem", str)
    @test occursin("LinearRightHandSide", str)
end

