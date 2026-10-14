# Executable companion to docs/src/guide.md.
# Run in an environment containing this package and the dependencies used below.

using NSDEBase
pendulum = IVP((u, t) -> [u[2], -sin(u[1])], [π/4, 0.0], (0.0, 10.0))
function lorenz!(du, u, t)
    du[1] = 10 * (u[2] - u[1])
    du[2] = u[1] * (28 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8/3) * u[3]
    return nothing
end
lorenz = IVP(lorenz!, [2.0, 3.0, -14.0], (0.0, 20.0))
@assert lorenz.rhs(lorenz.u0, 0.0) ≈ [10.0, 81.0, 130/3]
nothing

# u' = L*u: implicit solvers can use a direct linear solve.
L = [-1.0 0.0; 0.0 -100.0]
decay = IVP(L, [1.0, 1.0], (0, 1))

# u' = L*u + g(t): forcing takes time alone.
forced = IVP(LRHS(L, t -> [sin(t), 0.0]), [0.0, 0.0], (0, 1))

# u' = L*u + f_nonstiff(u,t): for IMEX or exponential methods.
split = IVP(L, (u, t) -> [sin(u[2]), cos(u[1])], [1.0, 0.0], (0, 1))
@assert decay.tspan === (0.0, 1.0)
@assert split.rhs isa SplitRightHandSide
nothing

integer_data = IVP((u, t) -> -u, [1, 2], (0, 1))
rotation = IVP(-1im, [1.0], (0, 1))
new_start = [1.0, 1.0, 1.0]
short = remake(lorenz; u0=new_start, tspan=(0.0, 2.0))
@assert integer_data.u0 isa Vector{Float64}
@assert eltype(rotation.u0) <: Complex
@assert short.rhs === lorenz.rhs
@assert short.u0 === new_start
nothing
