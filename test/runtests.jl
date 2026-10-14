using NSDEBase
using LinearAlgebra
using SparseArrays
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
function measure_remake_allocs(problem, u0, tspan)
    remake(problem; u0, tspan)
    return @allocated remake(problem; u0, tspan)
end

struct TestChunkSolver <: AbstractInitialValueSolver
    state::Vector{Float64}
end

struct UndocumentedIVP <: AbstractInitialValueProblem end

@testset "NSDEBase" begin
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
            du_ref = [0.3 * 0.7, -0.2 * 1.2]         # u(1-u) by hand
            @test rhs(u, t) ≈ du_ref
            du = similar(u)
            @test rhs(du, u, t) ≈ du_ref             # f! derived from f
            @test du ≈ du_ref
            J_ref = Diagonal(1.0 .- 2.0 .* u)        # d/du [u(1-u)] = 1 - 2u
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
            @test rhs(u, t) ≈ du_ref                 # f derived from f!
            J_ref = [u[2] u[1]; cos(u[1]) 0.0]       # Jacobian by hand
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
            # A mutating function owes its caller nothing but the filled buffer,
            # so the derived out-of-place `f` must return that buffer, not the
            # return value of `g!`.
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
            @test rhs.Df(u, t) ≈ λ * Matrix(I, 2, 2) atol=1e-7
            # Holomorphic case: f(u) = u², ∂f/∂u = 2u on the diagonal
            rhs2 = RHS((u, t) -> u .^ 2; iscomplex=true)
            @test rhs2.Df(u, t) ≈ Diagonal(2.0 .* u) atol=1e-7
        end

        @testset "Wirtinger block extraction" begin
            A = ComplexF64[1+2im 3-im; -2im 4+im]
            B = ComplexF64[2-im 1+im; 3+2im -1im]
            u = ComplexF64[0.4 + 0.7im, -0.2 + 0.9im]
            f(u, t) = A * u + B * conj.(u)
            @test NSDEBase._wirtinger_jacobian(f, u, 0.0) ≈ A atol=1e-9
            @test NSDEBase._wirtinger_jacobian((u, t) -> u .^ 2, u, 0.0) ≈ Diagonal(2u) atol=1e-9
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
            g!(dg, t) = (dg[1]=sin(t); dg[2]=cos(t); dg)  # allocation-free by hand
            du_ref = L * u + g(t)
            for rhs in (LRHS(L, g), LRHS(L, g!))              # g and g! entry forms
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

        @testset "derived in-place forcing returns dense filled buffers" begin
            rhs = LRHS(hcat(-1.0), (du, t) -> (du.=2.0; nothing))
            @test rhs([1.0], 0.0) == [1.0]
            @test rhs.g(0.0) == [2.0]
            complexrhs = LRHS(hcat(-1.0), (du, t) -> (du.=1im; nothing))
            @test complexrhs([1.0 + 0im], 0.0) == [-1.0 + 1im]
            sparse_rhs = LRHS(SparseArrays.sparse(hcat(-1.0)), (du, t) -> (du.=2.0; nothing))
            @test sparse_rhs.g(0.0) isa Vector{Float64}
            @test sparse_rhs([1.0], 0.0) == [1.0]
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
                SRHS(RHS((u, t) -> L * u + g(t)), fₙₛ),      # nonlinear stiff part
                SRHS((u, t) -> L * u + g(t), RHS(fₙₛ)),      # Function stiff part
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
            # An LRHS-with-g stiff part: the scratch `v` serves both parts
            stiff = LRHS(L, (dg, t) -> (dg[1]=t; dg[2]=2t; dg))
            rhs = SRHS(stiff, RHS((du, u, t) -> du .= sin.(u)))
            du, v = similar(u), similar(u)
            @test rhs(du, v, u, t) ≈ du_ref
            @test measure_rhs_allocs(rhs, du, v, u, t) == 0
        end
    end

    @testset "calling conventions" begin
        @testset "default arguments and explicit in-place choice" begin
            f(u, t, σ=10) = -σ .* u
            @test RHS(f)([1.0], 0.0) == [-10.0]
            g(t, σ=10) = [σ + t]
            @test LRHS(-1.0, g)([1.0], 0.0) == [9.0]
            mutating(du, u, t=0.0) = (du.=(-u); nothing)
            rhs = RHS(mutating; inplace=true)
            @test rhs([1.0], 0.0) == [-1.0]
            @test IVP(mutating, [1.0], (0.0, 1.0); inplace=true).rhs([1.0], 0.0) == [-1.0]
            typed(u::AbstractVector, t::Real) = -u
            @test RHS(typed)([1.0], 0.0) == [-1.0]
            @test_throws ArgumentError RHS(f; inplace=:yes)
            variadic(u, t, args...) = -u
            @test_throws ArgumentError RHS(variadic)
            @test RHS(variadic; inplace=false)([1.0], 0.0) == [-1.0]
        end

        @testset "validation and forwarding" begin
            out(u, t) = -u
            mut!(du, u, t) = (du.=(-u); nothing)
            @test_throws ArgumentError RHS(out; inplace=true)
            @test_throws ArgumentError RHS(mut!; inplace=false)
            @test_throws ArgumentError LRHS(-1.0, t -> [1.0]; inplace=true)
            @test_throws ArgumentError LRHS(-1.0, (v, t) -> fill!(v, 1); inplace=false)
            err = try
                RHS(u -> -u)
                nothing
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin("f(u, t)", sprint(showerror, err))
            @test !occursin("variadic", sprint(showerror, err))
            ambiguous!(du, u, t=0.0) = (du.=(-u); nothing)
            for u in (1.0, [1.0])
                for make in (
                    (f, u) -> IVP(f, u, (0, 1); inplace=true),
                    (f, u) -> IVP(f, u, 0, 1; inplace=true),
                )
                    p = make(ambiguous!, u)
                    @test p.rhs(p.u0, 0.0) == [-1.0]
                end
                for make in (
                    u -> IVP(out, ambiguous!, u, (0, 1); inplace=(false, true)),
                    u -> IVP(out, ambiguous!, u, 0, 1; inplace=(false, true)),
                    u -> IVP(-1.0, ambiguous!, u, 0, 1; inplace=true),
                )
                    p = make(u)
                    @test p.rhs(p.u0, 0.0) == [-2.0]
                end
            end
            @test SRHS(out, ambiguous!; inplace=(false, true))([1.0], 0.0) == [-2.0]
            @test SRHS(ambiguous!, RHS(out); inplace=(true, nothing))([1.0], 0.0) == [-2.0]
            @test_throws ArgumentError IVP(out, out, [1.0], (0, 1); inplace=(false, true))
            @test_throws ArgumentError IVP(out, out, [1.0], (0, 1); inplace=(true, false, true))
            @test RHS((u, t, args...) -> -u; inplace=false)([1.0], 0.0) == [-1.0]
            @test_throws ArgumentError RHS((a, b, c, d, args...) -> a; inplace=false)
        end
    end

    @testset "InitialValueProblem" begin
        @testset "construction forms" begin
            f!(du, u, t) = du .= u
            problem = IVP(f!, 0.5, 0.0, 1.0)            # scalar u0, (t0, tN)
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

        @testset "numeric types are normalised without copying floating states" begin
            for state in ([1, 2], [1//2, 2//3], Complex{Int}[1 + im, 2 - im])
                p = IVP((u, t) -> -u, state, (0, 1))
                @test p.u0 == float.(state)
                @test eltype(p.u0) == float(eltype(state))
                @test p.u0 !== state
                @test p.tspan === (0.0, 1.0)
            end
            for T in (Float32, Float64, BigFloat, ComplexF64)
                u = T[1, 2]
                p = IVP(-1.0, u, (0, 1.0))
                @test p.u0 === u
                @test remake(p; u0=u).u0 === u
            end
            @test IVP(-1.0, Float32[1], (0.0f0, 1.0f0)).tspan === (0.0f0, 1.0f0)
            @test IVP((u, t) -> -u, [1.0], 0, 1.0).tspan === (0.0, 1.0)
            @test IVP(-1.0, 1, (0, 1)).u0 == [1.0]
        end

        @testset "complex operators promote states; function Jacobians infer complexity" begin
            u = [1.0]
            for rhs in (-1im, hcat(-1im), LRHS(-1im), SRHS(-1im, (u, t) -> zero.(u)))
                p = IVP(rhs, u, (0, 1))
                @test p.u0 == ComplexF64[1]
                @test p.u0 !== u
            end
            @test Dahlquist(0.5; λ=-1 + 2im).u0 == [0.5 + 0im]
            @test IVP(-1im, (u, t) -> zero.(u), u, (0, 1)).u0 == ComplexF64[1]
            q = IVP((u, t) -> -1im .* u, ComplexF64[1], (0, 1))
            @test q.rhs.Df(q.u0, 0.0) ≈ hcat(-1im) atol=1e-10
            split = IVP((u, t) -> -1im .* u, (u, t) -> zero.(u), ComplexF64[1], (0, 1))
            @test split.rhs.fₛ.Df(split.u0, 0.0) ≈ hcat(-1im) atol=1e-10
        end

        @testset "remake semantics" begin
            problem = Logistic([0.3], (0.0, 2.0))
            u0 = [0.9]
            sub = remake(problem; u0, tspan=(0.5, 1.5))
            @test sub.rhs === problem.rhs               # rhs shared by identity, not rebuilt
            @test sub.u0 === u0
            @test sub.tspan == (0.5, 1.5)
            @test remake(problem; tspan=(0.5, 1.5)).u0 === problem.u0
            @test (@inferred remake(problem; u0, tspan=(0.5, 1.5))) isa typeof(problem)
            # Hot seam: building a chunk problem must stay cheap (one small struct)
            @test measure_remake_allocs(problem, u0, (0.5, 1.5)) ≤ 128
        end

        @testset "remake defaults" begin
            p = Logistic(0.3, (0.0, 1.0))
            q = remake(p)
            @test q.rhs === p.rhs
            @test q.u0 === p.u0
            @test q.tspan == p.tspan
            @test remake(p; u0=0.7).u0 == [0.7]
            @test remake(p; u0=0.7).tspan == p.tspan
        end

        @testset "mixed time endpoints in remake" begin
            p = IVP((u, t) -> u, [1.0], (0, 10))
            q = @inferred remake(p; tspan=(0, 2.5))
            @test q.tspan === (0.0, 2.5)
            @test q.rhs === p.rhs
            @test q.u0 === p.u0
            @test remake(p; tspan=(0, big"2.5")).tspan isa Tuple{BigFloat,BigFloat}
        end

        @testset "custom IVP without remake" begin
            @test_throws ArgumentError remake(UndocumentedIVP(); u0=[1.0], tspan=(0, 1))
        end
    end

    @testset "chunk_solver" begin
        solver = TestChunkSolver([1.0])
        chunk = chunk_solver(solver, 0.5)
        @test chunk.state == solver.state
        @test chunk.state !== solver.state
        chunk.state[1] = 2
        @test solver.state == [1.0]
    end

    @testset "ODEs" begin
        s = sin(π / 4)

        @testset "hand-checked derivatives at u0, t0" begin
            checks = (
                (Dahlquist(), [0.5]),                                # λ·u0
                (Logistic(), [0.25]),                                # u0(1−u0)
                (SimplePendulum(), [0.0, -s]),                       # [ω, −sin θ]
                (DoublePendulum(), [0.0, 0.0, -s, 0.0]),             # M\v with θ₁=θ₂
                (VanDerPol(), [0.0, -1.0]),                          # [0, −u₁]
                (Rössler(), [0.0, 2.0, 0.2]),                        # by substitution
                (Lorenz(), [10.0, 81.0, 6.0 + 112.0 / 3.0]),         # by substitution
            )
            for (problem, du_ref) in checks
                du = similar(problem.u0)
                @test problem.rhs(du, problem.u0, problem.tspan[1]) ≈ du_ref
            end
        end

        @testset "Lorenz96" begin
            problem = Lorenz96()                        # u0 = [1.01; ones(39)], F = 8
            du = similar(problem.u0)
            problem.rhs(du, problem.u0, 0.0)
            @test du[1] ≈ 6.99                          # (u₂−u₃₉)u₄₀ − u₁ + F
            @test du[2] ≈ 7.0                           # (u₃−u₄₀)u₁ − u₂ + F
            @test du[3] ≈ 6.99                          # (u₄−u₁)u₂ − u₃ + F
            @test du[40] ≈ 7.01                         # (u₁−u₃₈)u₃₉ − u₄₀ + F
            @test_throws ArgumentError Lorenz96(ones(3))
        end

        @testset "three-positional forms have defaults" begin
            @test Dahlquist(0.5, 0.0, 1.0) isa InitialValueProblem
            @test Logistic(0.5, 0.0, 1.0) isa InitialValueProblem
            @test SimplePendulum([π/4, 0.0], 0.0, 1.0) isa InitialValueProblem
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
end
